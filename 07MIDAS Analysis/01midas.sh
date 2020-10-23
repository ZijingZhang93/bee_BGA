#! /usr/bin/bash -w 
###########This part will remove reads belonging to host and trim reads.###########
#Prepare workspace
cd genetic_workspace
mkdir 10micro

#This step index Amel_HAv3.1 genomes to prepare for the KneadData database.
bowtie2-build AAscripts_ref/Amel_HAv3.1_genomic/Amel_HAv3.1_genomic.fa AAscripts_ref/Amel_HAv3.1_genomic/Amel_HAv3.1_genomic.fa
cd 10micro
for i in `ls ../01cleanfq/*R1_fastpout.fq.gz`
do
	j=`basename $i`
	#This step will perform principled in silico separation of bacterial reads from these “contaminant” reads, be they from the honeybee.
	kneaddata --input ${i} \
	--input ${i%R1_fastpout.fq.gz}R2_fastpout.fq.gz \
	--reference-db ./AAscripts_ref/Amel_HAv3.1_genomic/ \
	--output knea \
	--threads 16
	
	#This step will trim low-quality bases from both ends using the Phred algorithm.
	seqtk trimfq knea/${j}_kneaddata_paired_1.fastq >${j%R1_fastpout.fq.gz}_micro_trim_R1.fq
	seqtk trimfq knea/${j}_kneaddata_paired_2.fastq >${j%R1_fastpout.fq.gz}_micro_trim_R2.fq
done
############This part will quantify bacterial species abundance############
for i in $(cat sample.list)
do
	#Estimating the read depth and relative abundance of bacterial species
	run_midas.py species midas_result/${i} -1 ${i}_micro_trim_R1.fastq.gz -2 ${i}_micro_trim_R2.fastq.gz -d /home/wujiaqiang/midas_database/ -t 16
done
#This part will merge results across samples.
merge_midas.py species merge_species -i midas_result -t dir