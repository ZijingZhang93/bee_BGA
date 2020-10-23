import sys
import os

#Usage:Align the reads to the reference genome.
#Example:python 02hisat.py /home/muxiaohuan/RNA_analysis/test/data/sample

#save sample files in dictionary format
os.chdir(sys.argv[1])
sample = os.listdir(sys.argv[1])
path = {}
for seq in sample:
    name = seq.split('.')[0][:-2]
    if name not in path:
        path[name]=[]
    path[name].append(os.path.abspath(seq))

#get reference genome index
os.chdir('../ref')
ref = os.listdir('./')
prefix = ref[2][:-4]
index = os.path.abspath('../index/'+prefix)

os.chdir('../../')
if not(os.path.isdir('Hisat')):
    os.mkdir('Hisat')
os.chdir('./Hisat')

#mapping reads using hisat2
for k in path:
    read1 = path[k][0]
    read2 = path[k][1]
    bam = k+'.sort.bam'
    os.system('hisat2 -p 16 --dta -x '+index+' -1 '+read1+' -2 '+read2+' | '\
              +'samtools view -bS -f 0x2 | samtools sort -@ 5 -O bam -o '+bam)
    os.system('samtools index '+bam)
    
