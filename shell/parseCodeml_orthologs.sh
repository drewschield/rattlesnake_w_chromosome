echo 'chrom\tstart\tcrotalus\tdeinagkistrodon\tdnds\tdn\tds'
for result in ./codeml/*.txt; do
	name=$result
	iname=${name%.txt*} # removes .txt from name
	iname=${iname##*/}  # removes directory name
	chrom=`echo $iname | cut -d'_' -f1`
	start=`echo $iname | cut -d'_' -f2`
	dacu=${iname##*mRNA-1_}
	crot=${iname%_Dacu*}
	dnds=`python parse_codeml_output.py $result | tail -n+2 | awk '{print $3}'`
	dn=`python parse_codeml_output.py $result | tail -n+2 | awk '{print $4}'`
	ds=`python parse_codeml_output.py $result | tail -n+2 | awk '{print $5}'`
	echo "$chrom\t$start\t$crot\t$dacu\t$dnds\t$dn\t$ds"
done
