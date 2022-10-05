"""
scaffold_list_extractor.py

This script is based on the simple 'scaffold_extractor.py' script for pulling out specific sequences from large fasta files,
but first reads in a list of sequence names saved in a file (single column), and returns fasta with desired sequences.
This of course will only return sequences if the IDs match perfectly, so it is best to view sequence names prior to use. 

Usage: python scaffold_list_extractor.py <in.fasta> <seq.list> <out.fasta>
"""

import sys
from Bio import SeqIO

fasta = SeqIO.parse(sys.argv[1], "fasta")
outfasta = open(sys.argv[3], 'w')

chromlist = []

for line in open(sys.argv[2], 'r'):
	chrom = line.split()[0]
# 	print chrom
	chromlist.append(chrom)

for record in fasta:
	scaffold = record.id
# 	print scaffold
	sequence = record.seq
	for chrom in chromlist:
		if scaffold == chrom:
			outfasta.write('>'+str(scaffold)+'\n'+str(sequence)+'\n')