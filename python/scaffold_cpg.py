"""
usage: python scaffold_cpg.py <input.fasta> <output.txt>

"""

import sys
from Bio import SeqIO

out = open(sys.argv[2], 'w')
out.write('chrom'+'\t'+'CpG_content'+'\n')

for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
	id = seq_record.id
	seq = seq_record.seq
	stringseq = str(seq)
	total = len(stringseq)
	cg = stringseq.count('GC')
	n = stringseq.count('N')
	new_total = total - n
	cpg_content = cg/float(new_total)
	
	out.write(id+'\t'+str(cpg_content)+'\n')