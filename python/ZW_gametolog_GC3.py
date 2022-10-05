import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import combinations

# specify input parameters
input = sys.argv[1]
chrom = 'scaffold-Z'
start = input.split('scaffold-Z_')[1]
start = start.split('_')[0]
#zlog = input.split('crovir-transcript-')[1]
#zlog = zlog.split('.')[0]
#zlog = 'crovir-transcript-'+zlog
#wlog = input.split('scaffold-W')[1]
#wlog = wlog.split('_crovir')[0]
#wlog = 'scaffold-W'+wlog

# deal with internal, spurious stop codons from transcripts
codon_stop_array = ["TAG", "TGA", "TAA"]

zw = SeqIO.parse(sys.argv[1],'fasta')
for record1, record2 in combinations(zw, 2):
	# remove sequence at end so that nucleotides are multiples of 3
	sequence1 = str(record1.seq)
	sequence2 = str(record2.seq)
	if '.33' in str(float(len(sequence1))/float(3)):
		sequence1 = sequence1[:-1]
	if '.33' in str(float(len(sequence2))/float(3)):
		sequence2 = sequence2[:-1]
	if '.66' in str(float(len(sequence1))/float(3)):
		sequence1 = sequence1[:-2]
	if '.66' in str(float(len(sequence2))/float(3)):
		sequence2 = sequence2[:-2]
	# remove internal stop codons present in any position in one or both sequences from both sequences!
	tempRecordSeq1 = list(sequence1)
	tempRecordSeq2 = list(sequence2)
	count = 0
	for index in range(0, len(sequence1), 3):
		codon1 = sequence1[index:index+3]
		codon2 = sequence2[index:index+3]
		if codon1 in codon_stop_array or codon2 in codon_stop_array:
			# this is important, keeps track of how many codons have been removed!
			index = index - (3*count)
			del tempRecordSeq1[index:index+3]
			del tempRecordSeq2[index:index+3]
			count = count + 1
	sequence1 = Seq("".join(tempRecordSeq1))
	sequence2 = Seq("".join(tempRecordSeq2))

# calculate GC3 for Z and W gametologs
zgc3 = 0
zpos3 = 0
for index in range(0, len(sequence1),3):
	pos3 = sequence1[index+2]
	zpos3 = zpos3 + 1
	if str(pos3) == 'G':
		zgc3 = zgc3 + 1
	if str(pos3) == 'C':
		zgc3 = zgc3 + 1

zgc3_stat = float(zgc3)/float(zpos3)
zgc3_prop = "%.3f" % zgc3_stat

wgc3 = 0
wpos3 = 0
for index in range(0, len(sequence2),3):
	pos3 = sequence2[index+2]
	wpos3 = wpos3 + 1
	if str(pos3) == 'G':
		wgc3 = wgc3 + 1
	if str(pos3) == 'C':
		wgc3 = wgc3 + 1

wgc3_stat = float(wgc3)/float(wpos3)
wgc3_prop = "%.3f" % wgc3_stat

print "%s\t%s\t%s\t%s" % (chrom,start,zgc3_prop,wgc3_prop)
