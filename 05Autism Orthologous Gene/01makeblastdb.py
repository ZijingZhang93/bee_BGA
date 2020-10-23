#*-* coding:utf8 *-*
#usage:build the database using the sequence we input
#Example:python 01makeblastdb.py <the sequence that we want to build database>
import sys
import os

makeblastdb = '/home/share/biosoft/blast/makeblastdb' #command
fa = sys.argv[1]
dbtype = 'prot' #['nucl','prot']
out = 'bee_prot' #prefix
order = makeblastdb+' -in '+fa+' -dbtype '+dbtype+' -out '+out+' -parse_seqids'
print(order)
os.system(order)
