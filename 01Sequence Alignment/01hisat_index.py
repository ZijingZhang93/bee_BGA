import sys
import os

#Usage:Create an index of the genome according to the input reference genome path,
#      The input path include the reference genome files in gff and fasta format.
#      Extract the cut site information and exon information.
#Example:python 01hisat_index.py /home/muxiaohuan/RNA_analysis/test/data/ref 

os.chdir(sys.argv[1])
ref = os.listdir(sys.argv[1])

#Maybe adjust according to the your situation
gene_gtf = os.path.abspath(ref[2]) #Reference genome in gff format
genome_fa = os.path.abspath(ref[1]) #Reference genome in fasta format
prefix = ref[2][:-4]

#Extract the cut site information and save it to a file 'chrX.ss'
os.system('extract_splice_sites.py '+gene_gtf+'>'+prefix+'.ss')
#Extract exon information from gtf file and save it to file 'chrX.exon'
os.system('extract_exons.py '+gene_gtf+'>'+prefix+'.exon')

ss = os.path.abspath(prefix+'.ss')
exon = os.path.abspath(prefix+'.exon')

os.chdir('../')
if not(os.path.isdir('index')):
    os.mkdir('index')
os.chdir('./index')
#creat index using hisat2-build
os.system('hisat2-build'+' --ss '+ss+' --exon '+exon+' '+genome_fa+' '+prefix)

