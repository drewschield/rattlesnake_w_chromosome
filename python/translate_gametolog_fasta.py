import sys
from Bio import SeqIO
from Bio import Seq

fastaname = str(sys.argv[1])
fastaname = fastaname.split('/')[1]
fastaname = fastaname.split('.fna')[0]
outfasta = open('ZW_gametolog_seq/'+fastaname+'.faa','w')

for record in SeqIO.parse(sys.argv[1],'fasta'):
	aa = record.seq.translate()
	outfasta.write('>'+str(record.id)+'\n'+str(aa)+'\n')
