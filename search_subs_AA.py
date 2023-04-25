#Produces tab-delimited table of aminoacid substitutions from an alignment
#Usage: python3 search_subs.py aa-alignment.fasta

import sys
import os
from Bio import SeqIO

if __name__ == "__main__":
	#print("Start")
	seqs = []
	for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
		seqs.append(seq_record)
	head = "\t"
	for k in range(len(seqs)):
		head = head + str(seqs[k].id) + "\t"
	print(head + "\n")
	
	for i in range(len(str(seqs[0].seq))):
		flag = 0
		for seq in seqs:
			if flag == 1:
				continue
			if str(seq.seq).upper()[i] != str(seqs[0].seq).upper()[i] and str(seq.seq).upper()[i] != "X" and str(seqs[0].seq).upper()[i] != "X":
				subs = ""
				for j in range(len(seqs)):
					subs = subs + str(seqs[j].seq).upper()[i] + "\t"
				print(str(i+1)+"\t"+subs+"\n")
				flag = 1

