#Calculates melting temperature for short oligos using a simple formula
#Takes tab-delimited list of oligos (Oligo_name\tOligo_seq)
#Usage python3 Rodney.py oligo_list.tsv

import sys

with open(sys.argv[1], 'r') as infile:
	for line in infile.readlines():
		A = 0
		C = 0
		G = 0
		T = 0
		for el in line.split('\t')[1].strip('\n'):
			if el in ['A', 'a']:
				A +=1
			elif el in ['C', 'c']:
				C +=1
			elif el in ['G', 'g']:
				G +=1
			elif el in ['T', 't']:
				T += 1
		Tm = 2*A+2*T+4*C+4*G
		print(line.split('\t')[0], line.split('\t')[1].strip('\n'), Tm)