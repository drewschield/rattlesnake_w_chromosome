dir=$1
out=$2
for seq in $dir/*.faa; do
	name=$seq
	iname=${name%.faa*}		# removes .faa from name
	iname=${iname##*/}			# removes directory name
	echo performing clustal-omega alignment on $seq
	clustalo -i $seq -o $out/$iname.aln.faa
done 

