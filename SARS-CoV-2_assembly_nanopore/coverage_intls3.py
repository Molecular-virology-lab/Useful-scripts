import sys
import os
import statistics

with open(sys.argv[1], 'r') as intervals:
	intls_arr = intervals.readlines()
out_str = ''
cnt = 1
out_str += sys.argv[2]
out_str += '\t'
with open(sys.argv[2]+'_cov4.txt', 'r') as dpth:
	dpth_arr = dpth.readlines()
cov_arr = []
for dp in dpth_arr:
	cov_arr.append(int(dp.split()[2]))
if len(cov_arr) > 0:
	tot_med = statistics.median(cov_arr)
else:
	tot_med = 0
out_str += str(tot_med)
out_str += '\t'
for intl in intls_arr:
	int_left = int(intl.split()[1]) - 1
	int_right = int(intl.split()[2])
	tmp_int_arr = cov_arr[int_left:int_right]
	if len(tmp_int_arr) > 0:
		int_med = statistics.median(tmp_int_arr)
	else:
		int_med = 0
	if tot_med != 0:
		out_str += str(round(int_med/tot_med, 2))
	else:
		out_str += '0'
	out_str += '\t'
with open(sys.argv[2]+'_cov_meds.txt', 'w') as ouf:
	ouf.write(out_str + "\n")
print(out_str)
