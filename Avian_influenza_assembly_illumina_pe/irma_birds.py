import os
import subprocess



samples = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48']
#samples = ['01']

for sample in samples:
	os.system("~/soft/flu-amd/IRMA FLU-utr "+sample+"* "+sample)
	#reply = subprocess.run(['/home/lmv/soft/flu-amd/IRMA', 'FLU-utr', sample+'*', sample])
	os.system("cat "+sample+"/"+"*.fasta > "+sample+"/"+sample+".fas")
	os.system("bwa-mem2 index "+sample+"/"+sample+".fas")
	os.system("cat "+sample+"* "+" > "+sample+".fastq.gz")
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