"""
usage: python scaffold_gc.py <input.fasta> <output.txt>

"""

import sys
from Bio import SeqIO

out = open(sys.argv[2], 'w')
out.write('chrom'+'\t'+'GC_content'+'\n')

for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
	id = seq_record.id
	seq = seq_record.seq
	stringseq = str(seq)
	total = len(stringseq)
	c = stringseq.count('C')
	g = stringseq.count('G')
	n = stringseq.count('N')
	gc = g + c
	new_total = total - n
	gc_content = gc/float(new_total)
	
	out.write(id+'\t'+str(gc_content)+'\n')