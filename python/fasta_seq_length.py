import sys
from Bio import SeqIO

# total = 0

for i in SeqIO.parse(sys.argv[1], "fasta"):
	id = i.id
	id = id.strip()
	length = len(i.seq)
	info = id+'\t'+str(length)
	print info
# 	total = total + int(length)
# print "sum of scaffold lengths =", total
