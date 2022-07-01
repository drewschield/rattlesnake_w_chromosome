echo 'chrom\tstart\tZW_gametologs\tdnds\tdn\tds'
for result in ./codeml/*.txt; do
	name=$result
	iname=${name%.txt*} # removes .txt from name
	iname=${iname##*/}  # removes directory name
	chrom=`echo $iname | cut -d'_' -f1`
	start=`echo $iname | cut -d'_' -f2`
	dnds=`python parse_codeml_output.py $result | tail -n+2 | awk '{print $3}'`
	dn=`python parse_codeml_output.py $result | tail -n+2 | awk '{print $4}'`
	ds=`python parse_codeml_output.py $result | tail -n+2 | awk '{print $5}'`
	echo "$chrom\t$start\t$iname\t$dnds\t$dn\t$ds"
done

