#BWAcycle tool
#Read mapping for influenza sequencing on Illumina
import sys
import os
import math
import shutil
import statistics
import operator
from datetime import datetime
import base64
import argparse
import zipfile
import matplotlib.pyplot as ppl
import matplotlib.ticker as ticker


pth = os.path.dirname(os.path.abspath(__file__))
bam_readcount_path = pth+os.sep+'bam-readcount'+os.sep


def get_cons2(bam, sample_number, gene, ref_file, step):
	#Получаем консенсусную последовательность и статистику по позициям
	print("Extracting base counts from BAM-file with Bam-readcount for sample "+sample_number+" and gene "+gene+". Further output from Bam-readcount:\n")
	bamrc_cmd = bam_readcount_path+'bam-readcount -w 0 -f '+ref_file+' '+bam+' > '+sample_number+'_'+gene+'_'+step+'_readcount.txt'
	os.system(bamrc_cmd)
	print("Bam-readcount finished\n")
	with open(sample_number+'_'+gene+'_'+step+'_readcount.txt', 'r') as f:
		rc_str_arr = f.readlines()
	rc = {}
	for el in rc_str_arr:
		rc_arr = el.split()
		rc[int(rc_arr[1])] = {"A":rc_arr[5].split(':'), "C":rc_arr[6].split(':'), "G":rc_arr[7].split(':'), "T":rc_arr[8].split(':'), "INDEL":[], "COV":int(rc_arr[3]), "REF":rc_arr[2]}
		if len(rc_arr) > 10:
			ind_arr = []
			for elt in rc_arr[10:]:
				ind_arr.append(elt.split(':'))
			rc[int(rc_arr[1])]["INDEL"] = ind_arr
	print("Calculating consensus for sample "+sample_number+" and gene "+gene+"\n")
	con = ""
	dels = [0,0]
	probl_file = open(sample_number+'_'+gene+'_problems'+step+'.txt', 'w')
	ks = list(rc.keys())
	ks.sort()
	l_cov_arr = []
	for i in ks:
		if rc[i]["COV"] < 10:
			l_cov_arr.append(i)
			probl_file.write("Coverage is lower than 10 in position "+str(i)+". Please review basecount data.\n")
		if dels[0] > 0 and dels[1] > rc[i]["COV"]:
			dels[0]-=1
		else:
			pos = rc[i]["REF"]
			if int(rc[i]["A"][1]) > int(rc[i]["C"][1]) and int(rc[i]["A"][1]) > int(rc[i]["G"][1]) and int(rc[i]["A"][1]) > int(rc[i]["T"][1]):
				pos = "A"
			if int(rc[i]["C"][1]) > int(rc[i]["A"][1]) and int(rc[i]["C"][1]) > int(rc[i]["G"][1]) and int(rc[i]["C"][1]) > int(rc[i]["T"][1]):
				pos = "C"
			if int(rc[i]["G"][1]) > int(rc[i]["C"][1]) and int(rc[i]["G"][1]) > int(rc[i]["A"][1]) and int(rc[i]["G"][1]) > int(rc[i]["T"][1]):
				pos = "G"
			if int(rc[i]["T"][1]) > int(rc[i]["C"][1]) and int(rc[i]["T"][1]) > int(rc[i]["G"][1]) and int(rc[i]["T"][1]) > int(rc[i]["A"][1]):
				pos = "T"
			if len(rc[i]["INDEL"]) > 0:
				indel_cov = [0, int(rc[i]["INDEL"][0][1])]
				for j in range(len(rc[i]["INDEL"])):
					if int(rc[i]["INDEL"][j][1]) > indel_cov[1]:
						indel_cov = [j, int(rc[i]["INDEL"][j][1])]
				if rc[i]["INDEL"][indel_cov[0]][0][0] == "+" and float(indel_cov[1])/float(rc[i]["COV"]) > 0.5:
					pos += rc[i]["INDEL"][indel_cov[0]][0][1:]
					print("Found more than 50 percent of insertions in position "+str(i)+". Please review basecount data.\n")
					probl_file.write("Found more than 50 percent of insertions in position "+str(i)+". Please review basecount data.\n")
				if rc[i]["INDEL"][indel_cov[0]][0][0] == "-" and float(indel_cov[1])/float(rc[i]["COV"]) > 0.5:
					dels = [len(rc[i]["INDEL"][indel_cov[0]][0][1:]), int(rc[i]["INDEL"][indel_cov[0]][1])]
					print("Found more than 50 percent of deletions in position "+str(i)+". Please review basecount data.\n")
					probl_file.write("Found more than 50 percent of deletions in position "+str(i)+". Please review basecount data.\n")
			con += pos
	print("Coverage is lower than 10 in "+str(len(l_cov_arr))+" following positions:\n"+str(l_cov_arr))
	if step != '':
		out_file = open(sample_number+'_'+gene+'_'+step+'.fasta', 'w')
		out_file.write('>'+sample_number+'_'+gene+'_intermediate_cons'+'\n'+con+'\n')
	else:
		out_file = open(sample_number+'_'+gene+'.fasta', 'w')
		out_file.write('>'+sample_number+'_'+gene+'\n'+con+'\n')
	out_file.close()
	probl_file.close()
	rc["CON"] = con
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return rc

def get_stats(cov_arr, sample_number, gene, sample_name):
	#Получение консенсусной последовательности, таблиц и графиков из массива частот встречаемости
	print("Calculating statistics for sample "+sample_number+" and gene "+gene+"\n")
	t = open(sample_number+'_'+gene+'.table', 'w')
	t.write("Position\tA\tC\tG\tT\tpA\tpC\tpG\tpT\tShannon\tCoverage\tDiversity\tConsensus\tA_SNP\tC_SNP\tG_SNP\tT_SNP\n")
	cons = cov_arr.pop("CON")
	Shan = []
	Cov = []
	Div = []
	bq_arr = []
	mq_arr = []
	sln = len(cov_arr)
	prob_count = 0
	ks = list(cov_arr.keys())
	ks.sort()
	for j in ks:
		Acount = int(cov_arr[j]["A"][1])
		Gcount = int(cov_arr[j]["G"][1])
		Ccount = int(cov_arr[j]["C"][1])
		Tcount = int(cov_arr[j]["T"][1])
		CovJ = cov_arr[j]["COV"]
		cn = {"A":Acount, "C":Ccount, "G":Gcount, "T":Tcount}
		pos = max(cn, key=cn.get)
		bq_arr.append(float(cov_arr[j][pos][3]))
		mq_arr.append(float(cov_arr[j][pos][2]))
#Расчет энтропии Шеннона
		PA = 0.0
		PG = 0.0
		PC = 0.0
		PT = 0.0
		if float(CovJ) != 0.0:
			PA = float(Acount)/float(CovJ)
			PG = float(Gcount)/float(CovJ)
			PC = float(Ccount)/float(CovJ)
			PT = float(Tcount)/float(CovJ)
		AENT = 0.0
		GENT = 0.0
		CENT = 0.0
		TENT = 0.0
		if PA != 0.0:
			AENT = -1*PA*math.log2(PA)
		if PG != 0.0:
			GENT = -1*PG*math.log2(PG)
		if PC != 0.0:
			CENT = -1*PC*math.log2(PC)
		if PT != 0.0:
			TENT = -1*PT*math.log2(PT)
		ShanJ = AENT + GENT + CENT + TENT
		DivJ = 0
		if CovJ > 10:
			DivJ = (Acount*Ccount + Acount*Gcount + Acount*Tcount + Ccount*Gcount + Ccount*Tcount + Gcount*Tcount)/((CovJ*CovJ-CovJ)/2)
#Вывод результатов в таблицу
		t.write(str(j)+"\t"+str(Acount)+"\t"+str(Ccount)+"\t"+str(Gcount)+"\t"+str(Tcount)+"\t"+str(PA)+"\t"+str(PC)+"\t"+str(PG)+"\t"+str(PT)+"\t"+str(ShanJ)+"\t"+str(CovJ)+"\t"+str(DivJ)+"\t"+pos)
#Поиск SNP
		if (PA > 0.05 and PA < 0.95):
			t.write("\t"+str(PA))
		else:
			t.write("\t")
		if (PC > 0.05 and PC < 0.95):
			t.write("\t"+str(PC))
		else:
			t.write("\t")
		if (PG > 0.05 and PG < 0.95):
			t.write("\t"+str(PG))
		else:
			t.write("\t")
		if (PT > 0.05 and PT < 0.95):
			t.write("\t"+str(PT))
		else:
			t.write("\t")
		t.write("\n")
		Shan.append(ShanJ)
		Cov.append(CovJ)
		Div.append(DivJ)
	t.close()
	print("Calculation finished\n")
#Рисование графика энтропии Шеннона
	print("Drawing Shannon graph\n")
	mSh = statistics.mean(Shan)
	vSh = statistics.stdev(Shan)
	v2Sh = vSh*2
	v3Sh = vSh*3
	v5Sh = vSh*5
	v10Sh = vSh*10
	mvector = []
	vvector = []
	v2vector = []
	v3vector = []
	v5vector = []
	v10vector = []
	qua_lim_arr = []
	low_qua_lim_arr = []
	cov_lim_arr = []
	low_cov_lim_arr = []
	wth = len(ks)/1000*6.4
	for i in range(0, len(ks)):
		mvector.append(mSh)
		vvector.append(mSh+vSh)
		v2vector.append(mSh+v2Sh)
		v3vector.append(mSh+v3Sh)
		v5vector.append(mSh+v5Sh)
		v10vector.append(mSh+v10Sh)
		qua_lim_arr.append(30)
		cov_lim_arr.append(100)
		low_cov_lim_arr.append(10)
		low_qua_lim_arr.append(20)
	ppl.figure(figsize=(60,10), dpi=200)
	ppl.plot(ks, Shan, 'b', ks, mvector, 'g--', ks, v3vector, 'r--', ks, v5vector, 'y--', ks, v10vector, 'k--')
	ppl.title(sample_name+" "+gene+" Shannon enthropy, plus 3, 5, 10 stdev")
	ppl.savefig(sample_number+'_'+gene+'_Shannon.svg')
	ppl.savefig(sample_number+'_'+gene+'_Shannon.png')
	ppl.clf()
	ppl.cla()
	print("Shannon graph ready\n")
#Рисование графика покрытия
	print("\nDrawing coverage graph\n")
	m = round(statistics.median(Cov), 2)
	fig, ax = ppl.subplots()
	ax.plot(ks, Cov, 'b', ks, cov_lim_arr, 'y--', ks, low_cov_lim_arr, 'r--')
	ppl.title(sample_name+" "+gene+",\nmedian coverage = "+str(m))
	if len(ks)//10 > 1000:
		ax.xaxis.set_major_locator(ticker.MultipleLocator(round(len(ks)//10, -3)))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(round(len(ks)//100, -2)))
	elif len(ks)//10 > 100 and len(ks) <= 1000:
		ax.xaxis.set_major_locator(ticker.MultipleLocator(round(len(ks)//10, -2)))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(round(len(ks)//100, -1)))
	elif len(ks)//10 > 10 and len(ks) <= 100:
		ax.xaxis.set_major_locator(ticker.MultipleLocator(round(len(ks)//10, -1)))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(round(len(ks)//100, 0)))
	ax.grid(which='major', color = 'k')
	ax.minorticks_on()
	#ax.grid(which='minor', color = 'grey', linestyle = ':')
	ax.set_yscale('log')
	ax.set_xlim(0, len(Cov))
	ax.set_ylim(1, 10000)
	y_major_ticks = [el for el in ax.get_yticks() if el >=1 and el <= 10000]
	ax.set_yticks(y_major_ticks)
	ax.set_yticklabels([str(int(el)) for el in y_major_ticks])
	fig.set_figwidth(wth)
	fig.set_figheight(5)
	ppl.savefig(sample_number+'_'+gene+'_Coverage.svg')
	ppl.savefig(sample_number+'_'+gene+'_Coverage.png')
	print("\nCoverage graph ready\n")
	fig.clf()
	ax.cla()
	ppl.close(fig)
	#Рисование диаграммы долей покрытия
	print("\nDrawing coverage diagram\n")
	count0 = 0
	count10 = 0
	count100 = 0
	count1000 = 0
	for el in Cov:
		if el <=10:
			count0 +=1
		elif el >10 and el <=100:
			count10 +=1
		elif el >100 and el <=1000:
			count100 +=1
		elif el >1000:
			count1000 +=1
	data_names = ['<=10', '10-100', '100-1000', '>1000']
	data_values = [count0, count10, count100, count1000]
	ppl.pie(data_values, labels = data_names, autopct='%1.2f%%', radius = 1.1)
	ppl.title('Coverage distribution for sample: \n'+sample_name)
	ppl.axis('equal')
	ppl.savefig(sample_number+'_'+gene+'_cov_distribution.svg')
	ppl.savefig(sample_number+'_'+gene+'_cov_distribution.png')
	ppl.clf()
	ppl.cla()
	ppl.close()
#Рисование графика качества прочтения
	print("Drawing read quality graph\n")
	bqm = round(statistics.median(bq_arr), 2)
	#ppl.figure(figsize=(8,6), dpi=80)
	ppl.plot(ks, bq_arr, 'b', ks, qua_lim_arr, 'r--', ks, low_qua_lim_arr, 'y--')
	ppl.title(sample_name+" "+gene+", median base quality = "+str(bqm))
	ppl.savefig(sample_number+'_'+gene+'_BaseQuality.svg')
	ppl.savefig(sample_number+'_'+gene+'_BaseQuality.png')
	print("Quality graph ready\n")
	ppl.clf()
	ppl.cla()
#Рисование графика качества выравнивания
	print("Drawing mapping quality graph\n")
	mqm = round(statistics.median(mq_arr), 2)
	ppl.plot(ks, mq_arr, 'b', ks, qua_lim_arr, 'r--')
	ppl.title(sample_name+" "+gene+", median mapping quality = "+str(mqm))
	ppl.savefig(sample_number+'_'+gene+'_MapQuality.svg')
	ppl.savefig(sample_number+'_'+gene+'_MapQuality.png')
	print("Quality graph ready\n")
	ppl.clf()
	ppl.cla()
#Рисование графика неоднородности
	print("Drawing sequence diversity graph\n")
	sdm = statistics.mean(Div)
	vd = statistics.stdev(Div)
	sdm_vector = []
	sdmv_vector = []
	sdmv2_vector = []
	sdmv3_vector = []
	for i in range(0, len(ks)):
		sdm_vector.append(sdm)
		sdmv_vector.append(sdm + vd)
		sdmv2_vector.append(sdm + vd*2)
		sdmv3_vector.append(sdm + vd *3)
	ppl.plot(ks, Div, 'b', ks, sdm_vector, 'g--', ks, sdmv_vector, 'r--', ks, sdmv2_vector, 'y--',ks, sdmv3_vector, 'k--')
	ppl.title(sample_name+" "+gene+", mean sequence diversity = "+str(sdm))
	ppl.savefig(sample_number+'_'+gene+'_SeqDiversity.svg')
	ppl.savefig(sample_number+'_'+gene+'_SeqDiversity.png')
	print("Sequence diversity ready\n")
	ppl.clf()
	ppl.cla()
	with open(sample_number+'_'+gene+'_problems.txt', 'r') as probl_file:
		probs = probl_file.readlines()
		prob_count = len(probs)
	ans = []
	ans.append(len(cons))
	ans.append(m)
	ans.append(prob_count)
	ans.append(bqm)
	ans.append(mqm)
	ans.append(sdm)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return ans

descriptionText = "Illumina read mapping for influenza"
parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-bam", dest="bam", required="true", help="BAM-file")
parser.add_argument("-sample_number", dest="sample_number", required="true", help="Sample number")
parser.add_argument("-ref", dest="ref", required="true", help="Path to reference file")
parser.add_argument("-gene", dest="gene", required="true", help="Gene")
args = parser.parse_args()

cnt1 = get_cons2(args.bam, args.sample_number, args.gene, args.ref, '')
print("Length of gene "+args.gene+" of sample "+args.sample_number+" after second run: "+str(len(cnt1)-1)+"\n")
cv = get_stats(cnt1, args.sample_number, args.gene, args.sample_number)
print(cv)