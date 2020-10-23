#*-* coding:utf8 *-*
#usage:find the protein sequence in the database we build using 'makeblastdb'
#Example:python 02blastp.py <query seq> <subject db>

import sys
import os

blastp = '/home/share/biosoft/blast/blastp' #command
query = sys.argv[1]
subject_db = sys.argv[2]
out_file = 'blast.out' #output filename
order = blastp+' -query '+query+' -task blastp -db '+subject_db+' -out '+out_file+' -outfmt 6'
print(order)
os.system(order)
