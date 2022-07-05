"""
scaffold_mdg4_repeat_content.py

usage: python scaffold_mdg4_repeat_content.py <scaff_list_with_lengths> <repeat_GFF> <output>
"""

import sys

out = open(sys.argv[3], 'w')

out.write('scaffold'+'\t'+'length'+'\t'+'repeat_length'+'\t'+'prop_repeat'+'\t'+'Gypsy_length'+'\t'+'prop_Gypsy'+'\n')

for line in open(sys.argv[1], 'r'):
	chrom = line.split()[0]
	length = int(line.split()[1])

	tot_repeat = 0
	tot_gypsy = 0
	for entry in open(sys.argv[2], 'r'):
		if not entry.lstrip().startswith('#'):
			scaff = entry.split()[0]
			start = int(entry.split()[3])
			end = int(entry.split()[4])
			bases = int(abs(end - start))
			feat = entry.split()[8]
			
			if scaff == chrom:
				tot_repeat = tot_repeat + bases
				
				if "Gypsy" in str(feat):
					tot_gypsy = tot_gypsy + bases
	
	prop_repeat = float(tot_repeat)/float(length)
	if prop_repeat > 1.0:
		prop_repeat = 1.0
	prop_gypsy = float(tot_gypsy)/float(length)
	if prop_gypsy > 1.0:
		prop_gypsy = 1.0

	print chrom, length, tot_repeat, prop_repeat, tot_gypsy, prop_gypsy
	out.write(str(chrom)+'\t'+str(length)+'\t'+str(tot_repeat)+'\t'+str(prop_repeat)+'\t'+str(tot_gypsy)+'\t'+str(prop_gypsy)+'\n')

