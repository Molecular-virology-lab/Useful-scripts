import sys
import matplotlib.pyplot as plt


with open(sys.argv[1], 'r') as inf:
	genes = []
	reads = []
	for line in inf.readlines():
		if line.split('\t')[0] != "*":
			genes.append(line.split('\t')[0])
			reads.append(int(line.split('\t')[2]))
		else:
			genes.append("Unmapped")
			reads.append(int(line.split('\t')[3]))

with open(sys.argv[2], 'r') as inf:
	pos1 = []
	cov1 = []
	name1 = sys.argv[2].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos1.append(int(line.split('\t')[1]))
		cov1.append(int(line.split('\t')[2]))

with open(sys.argv[3], 'r') as inf:
	pos2 = []
	cov2 = []
	name2 = sys.argv[3].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos2.append(int(line.split('\t')[1]))
		cov2.append(int(line.split('\t')[2]))

with open(sys.argv[4], 'r') as inf:
	pos3 = []
	cov3 = []
	name3 = sys.argv[4].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos3.append(int(line.split('\t')[1]))
		cov3.append(int(line.split('\t')[2]))

with open(sys.argv[5], 'r') as inf:
	pos4 = []
	cov4 = []
	name4 = sys.argv[5].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos4.append(int(line.split('\t')[1]))
		cov4.append(int(line.split('\t')[2]))

with open(sys.argv[6], 'r') as inf:
	pos5 = []
	cov5 = []
	name5 = sys.argv[6].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos5.append(int(line.split('\t')[1]))
		cov5.append(int(line.split('\t')[2]))

with open(sys.argv[7], 'r') as inf:
	pos6 = []
	cov6 = []
	name6 = sys.argv[7].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos6.append(int(line.split('\t')[1]))
		cov6.append(int(line.split('\t')[2]))

with open(sys.argv[8], 'r') as inf:
	pos7 = []
	cov7 = []
	name7 = sys.argv[8].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos7.append(int(line.split('\t')[1]))
		cov7.append(int(line.split('\t')[2]))

with open(sys.argv[9], 'r') as inf:
	pos8 = []
	cov8 = []
	name8 = sys.argv[9].split('.')[0].split('_')[1]
	for line in inf.readlines():
		pos8.append(int(line.split('\t')[1]))
		cov8.append(int(line.split('\t')[2]))


fig, axs = plt.subplots(3, 3, figsize=(30, 30))
axs[0, 0].pie(reads, labels=genes, autopct='%1.2f%%', shadow=True,explode=(0, 0, 0, 0, 0, 0, 0, 0, 0.1))
axs[0, 0].set_title('Total reads: '+str(sum(reads)))
axs[0, 1].semilogy(pos1, cov1)
axs[0, 1].set_title(name1)
axs[0, 1].set_ylim(1,100000)
axs[0, 1].grid(True)
axs[0, 2].semilogy(pos2, cov2)
axs[0, 2].set_title(name2)
axs[0, 2].set_ylim(1,100000)
axs[0, 2].grid(True)
axs[1, 0].semilogy(pos3, cov3)
axs[1, 0].set_title(name3)
axs[1, 0].set_ylim(1,100000)
axs[1, 0].grid(True)
axs[1, 1].semilogy(pos4, cov4)
axs[1, 1].set_title(name4)
axs[1, 1].set_ylim(1,100000)
axs[1, 1].grid(True)
axs[1, 2].semilogy(pos5, cov5)
axs[1, 2].set_title(name5)
axs[1, 2].set_ylim(1,100000)
axs[1, 2].grid(True)
axs[2, 0].semilogy(pos6, cov6)
axs[2, 0].set_title(name6)
axs[2, 0].set_ylim(1,100000)
axs[2, 0].grid(True)
axs[2, 1].semilogy(pos7, cov7)
axs[2, 1].set_title(name7)
axs[2, 1].set_ylim(1,100000)
axs[2, 1].grid(True)
axs[2, 2].semilogy(pos8, cov8)
axs[2, 2].set_title(name8)
axs[2, 2].set_ylim(1,100000)
axs[2, 2].grid(True)
outname = sys.argv[1].split('.')[0].split('_')[0]
plt.savefig(outname+'.png')