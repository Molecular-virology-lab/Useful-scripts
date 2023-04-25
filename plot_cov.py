#Plots coverage diagram from samtools depth output
#Usage:  samtools depth -aa -d0 file.bam > depth.txt && python3 plot_cov.py depth.txt
import sys
from pylab import *
import statistics
pos = []
cov = []
zeros = 0
with open(sys.argv[1], 'r') as inf:
	for line in inf:
		s = line.split()
		pos.append(int(s[1]))
		cov.append(int(s[2]))
		if int(s[2]) == 0:
			zeros +=1
		#m = m + int(s[2])
	#m = m/pos[end]
	semilogy(pos, cov)
	savefig(sys.argv[1].split('.')[0]+'.png')
	print("Median: ", statistics.median(cov))
	print("Zeros: ", zeros)
		