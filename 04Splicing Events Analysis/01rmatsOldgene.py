import sys
import os

#Usage:Analysis of event-level differential splicing
#      based on the bee reference genome.
#Example:python 01rmatsOldgene.py /home/muxiaohuan/RNA_analysis/test/Hisat
#        python 01rmatsOldgene.py <Bam file after alighment>

#save sample files in dictionary format
dic = {}
os.chdir(sys.argv[1])
bam = os.listdir(sys.argv[1])
dic = {}
for f in bam:
    if f.endswith('.bam'):
        name_all = f.split('.')[0]
        name = name_all.split('_')[0][:-1]
        if name not in dic:
            dic[name] = []
        dic[name].append(os.path.abspath(f))

#get reference genome index
os.chdir('../data/ref/')
ref = os.listdir('./')
ref = os.path.abspath(ref[2])

os.chdir('../../')
if not(os.path.isdir('old_Rmats')):
    os.mkdir('old_Rmats')
os.chdir('./old_Rmats')
for k in dic:
    f = open(k+'.txt','w+')
    s = ''
    for b in dic[k]:
        s = s+b+','
    s = s.rstrip(',')
    f.write(s)
    f.close()

control = 'GF'
control_path = os.path.abspath(control+'.txt')
test = list(dic.keys())
test.remove(control)
dic_rmats = {}
for i in test:
    vs = i+'-'+control
    if vs not in dic_rmats:
        dic_rmats[vs] = []
    dic_rmats[vs].append(os.path.abspath(i+'.txt'))

#analysis of splicing events using rmats.py
for k in dic_rmats:
    os.system('echo'+' '+k)
    os.system('python ~/biosoft/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --b1 %s --b2 %s --gtf %s --od %s -t paired --readLength  %i --cstat %f --nthread %i'%(dic_rmats[k][0],control_path,ref,k,150,0.0001,10))
