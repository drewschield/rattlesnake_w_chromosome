import sys
from Bio import SeqIO

"""
Target files are:

annotation table: gametologs.cds.one2one.Z_geneID_coordinates.txt
Z transcripts: chrZ.cds.fasta
W transcripts: chrW.cds.fasta
"""

for line in open(sys.argv[1],'r'):
	chrom = line.split()[0]
	start = line.split()[1]
	wtran = line.split()[3]
	ztran = line.split()[4]
	
	fastaname = str(chrom+'_'+start+'_'+wtran+'_'+ztran)
	outfasta = open('ZW_gametolog_seq/'+fastaname+'.fna','w')
	
	for zrecord in SeqIO.parse(sys.argv[2],'fasta'):
		if str(zrecord.id) == str(ztran):
			outfasta.write('>crotalusZgametolog'+'\n'+str(zrecord.seq)+'\n')
	for wrecord in SeqIO.parse(sys.argv[3],'fasta'):
		if str(wrecord.id) == str(wtran):
			outfasta.write('>crotalusWgametolog'+'\n'+str(wrecord.seq)+'\n')

