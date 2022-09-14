import os
import subprocess



#samples = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '80', '81', '82', '88', '89', '90', '91', '92', '93', '94', '95', '96']
#samples = ['01']
samples = []
with open('sample_list.csv', 'r') as smp_list:
	for line in smp_list.readlines():
		samples.append(line.strip('\n'))

for sample in samples:
	print(sample)
	os.system("cat "+sample+"_* "+"> "+sample+".fastq.gz")
	os.system("IRMA FLU-utr "+sample+".fastq.gz "+sample)
	#reply = subprocess.run(['/home/lmv/soft/flu-amd/IRMA', 'FLU-utr', sample+'*', sample])
	os.system("cat "+sample+"/"+"*.fasta > "+sample+"/"+sample+".fas")
	os.system("bwa-mem2 index "+sample+"/"+sample+".fas")
	#os.system("cat "+sample+"* "+" > "+sample+".fastq.gz")
	os.system("bwa-mem2 mem -Y -t 4 "+sample+"/"+sample+".fas "+sample+".fastq.gz | samtools view -b | samtools sort -o "+sample+"/"+sample+".bam")
	os.system("rm "+sample+".fastq.gz")
	os.system("samtools index "+sample+"/"+sample+".bam")
	os.system("samtools idxstats "+sample+"/"+sample+".bam")
	os.system("samtools idxstats "+sample+"/"+sample+".bam > "+sample+"/"+sample+"_idxstats.txt")
	os.system("echo "+sample+" >> idxstats.txt")
	os.system("samtools idxstats "+sample+"/"+sample+".bam >> idxstats.txt")
	#result = subprocess.run(['ls', sample+'/'+'*.fasta'], stdout=subprocess.PIPE, encoding='utf-8')
	genes = []
	with open(sample+"/"+sample+".fas", "r") as inf:
		for line in inf.readlines():
			if line[0] == ">":
				genes.append(line[1:len(line)-1])
	depths = ""
	for segment in genes:
	#for gene in result.stdout:
		#segment = gene.split(".")[0]
		os.system("echo "+segment)
		os.system("samtools depth -aa -r "+segment+" "+sample+"/"+sample+".bam > "+sample+"/"+segment+"_depth.txt")
		os.system("samtools mpileup -aa -A -d 0 -Q 0 -r "+segment+" "+sample+"/"+sample+".bam | ivar consensus -p "+sample+"_"+segment+" -q 25 -m 3 -i "+sample+"_"+segment)
		depths += sample+"/"+segment+"_depth.txt "
		os.system("rm "+sample+"_"+segment+".qual.txt")
	os.system("python3 influenza_cov.py "+sample+"/"+sample+"_idxstats.txt "+depths)
	os.system("cp "+sample+"/"+sample+".png .")