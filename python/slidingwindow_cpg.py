"""
usage: python slidingwindow_cpg.py <input.fasta> <windowsize> <windowstep> <output.txt>

example usage to look at CpG content in 10Kb, non-overlapping windows:

python slidingwindow_cpg.py input.fasta 10000 10000 output.10Kb.txt

"""

import sys
from Bio import SeqIO

out = open(sys.argv[4], 'w')
out.write('chrom'+'\t'+'start'+'\t'+'end'+'\t'+'CpG_content'+'\n')

def window(seq, win, step):
	seqlen = len(seq)
	for i in range(0,seqlen,step):
		j = seqlen if i+win>seqlen else i+win
		yield seq[i:j]
		if j==seqlen: break

for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
	start = 1
	end = int(sys.argv[2])
	id = seq_record.id
	seq = seq_record.seq
	for subseq in window(seq, int(sys.argv[2]),int(sys.argv[3])):
		stringseq = str(subseq)
		total = len(stringseq)
		cg = stringseq.count('CG')
		n = stringseq.count('N')
		new_total = total - n
		if new_total == 0:
			cpg_content = 'NA'
		else:
			cpg_content = cg/float(new_total)
		
		out.write(id+'\t'+str(start)+'\t'+str(end)+'\t'+str(cpg_content)+'\n')
		end += int(sys.argv[2])
		start += int(sys.argv[2])
