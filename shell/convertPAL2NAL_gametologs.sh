dir=$1
seqdir=$2
for aln in $dir/*.aln.faa; do
	name=$aln
	iname=${name%.aln.faa*}		# removes .aln.faa from name
	iname=${iname##*/}			# removes directory name
	echo performing codon-based nucleotide alignment on $aln
	pal2nal.pl $aln $seqdir/$iname.fna -output paml -nogap > $dir/$iname.pal2nal
done 

