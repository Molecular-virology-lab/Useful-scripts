import sys
import os
from Bio import SeqIO


if __name__ == "__main__":
	seqs = []
	os.system('cat wuhan3524.fas '+sys.argv[1]+'_minimap_con.fasta '+'artic-nanopolish-cpu'+os.sep+sys.argv[2]+'.consensus.fasta '+'artic-medaka'+os.sep+sys.argv[2]+'.consensus.fasta > '+sys.argv[1]+'_con_all.fasta')
	os.system('mafft --auto --anysymbol --nomemsave '+sys.argv[1]+'_con_all.fasta > '+sys.argv[1]+'_con_aln.fasta')
	for seq_record in SeqIO.parse(sys.argv[1]+'_con_aln.fasta', 'fasta'):
		seqs.append(seq_record)
#	for seq_record in SeqIO.parse(sys.argv[1]+os.sep+sys.argv[1]+'_minimap_con.fasta', 'fasta'):
#		seqs.append(seq_record)
#	for seq_record in SeqIO.parse(sys.argv[1]+os.sep+'artic-nanopolish-cpu'+os.sep+sys.argv[2]+'.consensus.fasta', 'fasta'):
#		seqs.append(seq_record)
#	for seq_record in SeqIO.parse(sys.argv[1]+os.sep+'artic-medaka'+os.sep+sys.argv[2]+'.consensus.fasta', 'fasta'):
#		seqs.append(seq_record)
	with open(sys.argv[1]+'_minimap_cov.txt', 'r') as covf:
		cov = covf.readlines()
		cov_hash = {}
		for j in range(1,len(cov)):
			cov_s = cov[j].split()
			cov_hash[int(cov_s[0])] = cov_s
	out_l = open(sys.argv[1]+"_ltr_comp.txt", "w")
	out_n = open(sys.argv[1]+"_Ns_comp.txt", "w")
	out_gap = open(sys.argv[1]+"_gaps_comp.txt", "w")
	out_diff = open(sys.argv[1]+"_diffs_comp.txt", "w")
	
	out_l.write("Position\tMN908947.3\t3524_GISAID\tMinimap2\tARTIC_nanopolish\tARTIC_medaka\tMinimap_cov\tA\tC\tG\tT\n")
	out_n.write("Position\tMN908947.3\t3524_GISAID\tMinimap2\tARTIC_nanopolish\tARTIC_medaka\tMinimap_cov\tA\tC\tG\tT\n")
	out_gap.write("Position\tMN908947.3\t3524_GISAID\tMinimap2\tARTIC_nanopolish\tARTIC_medaka\tMinimap_cov\tA\tC\tG\tT\n")
	out_diff.write("Position\tMN908947.3\t3524_GISAID\tMinimap2\tARTIC_nanopolish\tARTIC_medaka\tMinimap_cov\tA\tC\tG\tT\n")
	
	if len(sys.argv) > 3:
		out_f = open(sys.argv[1]+"_final_con.fasta", "w")
		out_f.write(">"+sys.argv[2]+"\n")
	for i in range(len(seqs[0])):
		if len(seqs[3]) == i:
			seqs[3].seq +="N"
		if len(seqs[4]) == i:
			seqs[4].seq +="N"
		if i+1 in cov_hash.keys():
			pst = str(i+1)+"\t"+seqs[0][i].upper()+"\t"+seqs[1][i].upper()+"\t"+seqs[2][i].upper()+"\t"+seqs[3][i].upper()+"\t"+seqs[4][i].upper()+"\t"+cov_hash[i+1][1]+"\t"+cov_hash[i+1][2]+"\t"+cov_hash[i+1][3]+"\t"+cov_hash[i+1][4]+"\t"+cov_hash[i+1][10]+"\n"
		else:
			pst = str(i+1)+"\t"+seqs[0][i].upper()+"\t"+seqs[1][i].upper()+"\t"+seqs[2][i].upper()+"\t"+seqs[3][i].upper()+"\t"+seqs[4][i].upper()+"\t0\t0\t0\t0\t0\n"
		if (seqs[0][i].upper() != seqs[2][i].upper() or seqs[0][i].upper() != seqs[3][i].upper() or seqs[0][i].upper() != seqs[4][i].upper() or seqs[2][i].upper() != seqs[3][i].upper() or seqs[2][i].upper() != seqs[4][i].upper() or seqs[3][i].upper() != seqs[4][i].upper()) and (seqs[0][i].upper() != "N" and seqs[1][i].upper() != "N" and seqs[2][i].upper() != "N" and seqs[3][i].upper() != "N" and seqs[4][i].upper() != "N" and seqs[0][i].upper() != "-" and seqs[1][i].upper() != "-" and seqs[2][i].upper() != "-" and seqs[3][i].upper() != "-" and seqs[4][i].upper() != "-"):
			out_l.write(pst)
		if ((seqs[3][i].upper() == "N" or seqs[4][i].upper() == "N") and seqs[2][i] != "-"):
			out_n.write(pst)
		if (seqs[3][i].upper() == "N" and seqs[4][i].upper() == "N" and seqs[2][i] == "-"):
			out_gap.write(pst)
		if (seqs[3][i].upper() != seqs[4][i].upper()):
			out_diff.write(pst)
		if (seqs[3][i].upper() == "N" and seqs[4][i].upper() == "N" and seqs[2][i] != "-" and seqs[2][i].upper() != seqs[0][i].upper()):
			out_diff.write(pst)
		if len(sys.argv) > 3:
			if seqs[3][i].upper() != "N":
				out_f.write(seqs[3][i].upper())
			elif seqs[3][i].upper() == "N" and seqs[4][i].upper() != "N":
				out_f.write(seqs[4][i].upper())
			elif seqs[3][i].upper() == "N" and seqs[4][i].upper() == "N" and seqs[2][i].upper() != "-" and seqs[2][i].upper() == seqs[0][i].upper():
				out_f.write(seqs[2][i].upper())
			elif seqs[3][i].upper() == "N" and seqs[4][i].upper() == "N" and seqs[2][i].upper() != "-" and seqs[2][i].upper() != seqs[0][i].upper() and int(cov_hash[i+1][10]) > 50:
				out_f.write(seqs[2][i].upper())
			else:
				out_f.write("N")
	out_l.close()
	out_n.close()
	out_gap.close()
	out_diff.close()
	if len(sys.argv) > 3:
		out_f.close()