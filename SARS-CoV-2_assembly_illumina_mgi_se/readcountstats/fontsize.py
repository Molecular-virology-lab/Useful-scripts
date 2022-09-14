import sys
import os
from Bio import SeqIO
import statistics
seqs = []
for seq_record in SeqIO.parse(sys.argv[1]+'_tanoti_con.fa', 'fasta'):
	seqs.append(seq_record)
with open(sys.argv[1]+'_tanoti_depth.txt', 'r') as dpth:
		dpth_lines = dpth.readlines()
seq = str(seqs[0].seq)
nstr = ''
dpth_hash = {}
dpth_arr = []
for line in dpth_lines:
	line_split = line.split('\t')
	dpth_hash[int(line.split('\t')[1])-1] = int(line_split[2])
	dpth_arr.append(int(line_split[2]))
for i in range(len(seq)):
	if i in dpth_hash.keys() and dpth_hash[i] >= 10:
		nstr += seq[i].upper()
	else:
		nstr += seq[i].lower()
med = statistics.median(dpth_arr)
Acount = nstr.count('A') + nstr.count('a')
Ccount = nstr.count('C') + nstr.count('c')
Gcount = nstr.count('G') + nstr.count('g')
Tcount = nstr.count('T') + nstr.count('t')
Ncount = nstr.count('N') + nstr.count('n')
with open(sys.argv[1]+'_con.fasta', 'w') as ouf:
	ouf.write('>'+seqs[0].id+'_med_'+str(med)+'_ACTG_'+str(Acount+Ccount+Gcount+Tcount)+'_N_'+str(Ncount)+'\n'+nstr+'\n')
print(seqs[0].id+'_med_'+str(med)+'_ACTG_'+str(Acount+Ccount+Gcount+Tcount)+'_N_'+str(Ncount))
