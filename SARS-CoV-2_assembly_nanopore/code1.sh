#!/bin/bash
inputpath=/home/lmv/sdb1/20220825-2
touch stats.txt
touch cov_meds.txt
for sample in 01 02 03 04 05 06 07 08 09 10 11 12;
#for sample in 01;
do 
	mkdir $sample
	#cd 'barcode'$sample
	cd $sample
	echo $sample
	#cp '../../FASTQ/'$sample'/'$sample'.fastq.gz' $sample'.fastq.gz'
	cat $inputpath/barcode$sample/* > $sample.fastq.gz
	minimap2 -a -x map-ont ../WIV04.fasta $sample'.fastq.gz' | samtools view -bS | samtools sort -o $sample'.bam'
	samtools index $sample'.bam'
	samtools idxstats $sample'.bam'
	samtools depth -aa -d0 $sample'.bam' > $sample'_cov1.txt'
	ivar trim -i $sample'.bam' -b ../V4.bed -p $sample'_trimmed' -m 50 -q 9 -s 50 -e
	samtools sort -o $sample'_trimmed_sorted.bam' $sample'_trimmed.bam'
	samtools index $sample'_trimmed_sorted.bam'
	samtools depth -aa -d0 $sample'_trimmed_sorted.bam' > $sample'_cov2.txt'
	samtools mpileup -aa -A -d 0 -Q 0 -r WIV04_402124_2019-12-30:1-29869 $sample'_trimmed_sorted.bam' | ivar consensus -p $sample -q 9 -m 5 -i $sample
	minimap2 -a -x map-ont $sample.fa $sample'.fastq.gz' | samtools view -bS | samtools sort -o $sample'_2con.bam'
	samtools index $sample'_2con.bam'
	samtools idxstats $sample'_2con.bam'
	samtools depth -aa -d0 $sample'_2con.bam' > $sample'_cov3.txt'
	ivar trim -i $sample'_2con.bam' -b ../V4.bed -p $sample'_2con_trimmed' -m 50 -q 9 -s 50 -e
	samtools sort -o $sample'_2con_trimmed_sorted.bam' $sample'_2con_trimmed.bam'
	samtools index $sample'_2con_trimmed_sorted.bam'
	samtools depth -aa -d0 $sample'_2con_trimmed_sorted.bam' > $sample'_cov4.txt'
	samtools mpileup -aa -A -d 0 -Q 0 -r $sample':1-29860' $sample'_2con_trimmed_sorted.bam' | ivar consensus -p $sample'_final' -q 9 -m 10 -i $sample
	#medaka consensus --model r941_min_sup_variant_g507 $sample'_2con_trimmed_sorted.bam' $sample'_medaka_cons'
	medaka consensus --model r941_min_sup_g507  --threads 3 $sample'_2con_trimmed_sorted.bam' $sample'_medaka_cons'
	medaka stitch --threads 3 $sample'_medaka_cons' $sample'_final.fa' $sample'_final.fasta'
	#medaka variant $sample'_final.fa' $sample'_medaka_cons' $sample'_out.vcf'
	#medaka tools annotate --dpsp $sample'_out.vcf' $sample'_final.fa' $sample'_2con_trimmed_sorted.bam' $sample'_annotated.vcf'
	#bgzip $sample'_annotated.vcf'
	#tabix -p vcf $sample'_annotated.vcf.gz'
	#bcftools consensus -f $sample'_final.fa' -o $sample'_final.fasta' $sample'_annotated.vcf.gz'
	python3 ../fontsize_ont.py $sample'_final.fasta' $sample'_cov4.txt' >> ../stats.txt
	minimap2 -aY -x map-ont ../check.fas $sample'.fastq.gz' | samtools view -bS | samtools sort -o $sample'_check.bam'
	samtools index $sample'_check.bam'
	samtools idxstats $sample'_check.bam' >> ../stats.txt
	rm $sample'_check.bam'
	rm $sample'_check.bam.bai'
	cp $sample'_final_con.fasta' ../
	#cp $sample'_final.fa' ../$sample'_final.fasta'
	python3 ../plot_depth1.py -i $sample'_2con_trimmed_sorted.bam' -o $sample'_ito' -p ../V4.bed -r $sample'.fa' -t 16
	samtools mpileup -d 0 -f $sample'.fa' $sample'_2con_trimmed_sorted.bam' | python3 ../plot_indels_logy.py -o $sample'_indels.png'
	python3 ../coverage_intls3.py ../V4_nonoverlaps.txt $sample >> ../cov_meds.txt
	cp $sample'_ito.png' ../
	cp $sample'_indels.png' ../
	cd ..
done;
