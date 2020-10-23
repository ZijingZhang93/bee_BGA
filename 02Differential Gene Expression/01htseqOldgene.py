import sys
import os

#Usage:Quantification of gene expression
#Example:python 01htseqOldgene.py /home/muxiaohuan/RNA_analysis/test/Hisat 

#save sample files in dictionary format
os.chdir(sys.argv[1])
bam = os.listdir(sys.argv[1])
sample = {}
for f in bam:
    if f.endswith('.bam'):
        name_all = f.split('.')[0]
        if name_all not in sample:
            sample[name_all] = []
        sample[name_all].append(os.path.abspath(f))
        sample[name_all].append(name_all+'_htcount.txt')

#get reference genome in gff format
os.chdir('../data/ref')
ref = os.listdir('./')
ref_gtf = os.path.abspath(ref[0])

os.chdir('../../')
if not(os.path.isdir('old_HisatHtseq')):
    os.mkdir('old_HisatHtseq')
os.chdir('./old_HisatHtseq')

#count the expression of genes for every sample using htseq
for k in sample:
    os.system('/usr/local/anaconda/bin/htseq-count -f bam -r name -a 10 -t gene -s no -m union -i Name %s %s >%s'%(sample[k][0],ref_gtf,sample[k][1]))

#merge the htseq results
os.system('merge.py *_htcount.txt>merge_htcount.txt')
