#!/usr/bin/python

#Usage = [RBH.py BLASTOUTPUT1 BLASTOUTPUT2 outfile] 
import sys, re


infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile = sys.argv[3]

#parse first BLAST results
FL1 = open(infile1, 'r')
D1 = {} 
for Line in FL1:
	if ( Line[0] != '#' ):
		Line.strip()
		Elements = re.split('\t', Line)
		queryId = Elements[0]
		subjectId = Elements[1]
		if ( not ( queryId in D1.keys() ) ):
			D1[queryId] = subjectId  

#parse second BLAST results
FL2 = open(infile2, 'r')
D2 = {}
for Line in FL2:
	if ( Line[0] != '#' ):
		Line.strip()
		Elements = re.split('\t', Line)
		queryId = Elements[0]
		subjectId = Elements[1]
		if ( not ( queryId in D2.keys() ) ):
			D2[queryId] = subjectId  

#Now, pick the share pairs

SharedPairs={}
for id1 in D1.keys():
	value1 = D1[id1]
	if ( value1 in D2.keys() ):
		if ( id1 == D2[value1] ) : 
			SharedPairs[value1] = id1


outfl = open( outfile, 'w')

for k1 in SharedPairs.keys():
	line = k1 + ',' + SharedPairs[k1] + '\n'
	outfl.write(line)

outfl.close()

print"Done. Reciprocal best Hits from", sys.argv[1], "and", sys.argv[2], "are in", sys.argv[3]
