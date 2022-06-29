import sys
from Bio import SeqIO

"""
Target files are:

annotation table: orthologs_crotalus_anolis.one2one.txt
anolis CDS: anolis.cds.fasta
crotalus CDS: crotalus.cds.fasta
"""

ortholog = 1
for line in open(sys.argv[1],'r'):
	crot = line.split(",")[0]
	anol = line.split(",")[1]
	anol = anol.strip()
	
	fastaname = str('crotalus'+'_'+'anolis'+'_ortholog'+str(ortholog))
	print fastaname
	outfasta = open('ortholog_sequences/'+fastaname+'.fna','w')
	
	for crotrecord in SeqIO.parse(sys.argv[2],'fasta'):
		if str(crotrecord.id) == str(crot):
			outfasta.write('>'+str(crotrecord.id)+'\n'+str(crotrecord.seq)+'\n')
	for anolrecord in SeqIO.parse(sys.argv[3],'fasta'):
		if str(anolrecord.id) == str(anol):
			outfasta.write('>'+str(anolrecord.id)+'\n'+str(anolrecord.seq)+'\n')
	ortholog = ortholog + 1

