import sys
from Bio import SeqIO

"""
Target files are:

annotation table: crotalusZ_deinagkistrodon.ZW_gametolog.coordinates.txt
Z CDS: chrZ.cds.fasta
Deinagkistrodon CDS: deinagkistrodon_acutus.cds.fasta
"""

for line in open(sys.argv[1],'r'):
	chrom = line.split()[0]
	start = line.split()[1]
	crot = line.split()[3]
	dein = line.split()[4]
	
	fastaname = str(chrom+'_'+start+'_'+crot+'_'+dein)
	outfasta = open('Z_ortholog_seq/'+fastaname+'.fna','w')
	
	for crecord in SeqIO.parse(sys.argv[2],'fasta'):
		if str(crecord.id) == str(crot):
			outfasta.write('>crotalusZ'+'\n'+str(crecord.seq)+'\n')
	for drecord in SeqIO.parse(sys.argv[3],'fasta'):
		if str(drecord.id) == str(dein):
			outfasta.write('>deinagkistrodonZ'+'\n'+str(drecord.seq)+'\n')
