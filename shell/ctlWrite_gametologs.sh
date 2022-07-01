dir=$1
ctldir=$2
codemldir=$3
for aln in $dir/*.pal2nal; do
	name=$aln
	iname=${name%.pal2nal*}		# removes pal2nal from name
	iname=${iname##*/}			# removes directory name
	scaff=`echo $iname | cut -d '_' -f 1`
	start=`echo $iname | cut -d '_' -f 2`
	abbrev=`echo ${scaff}_${start}`
	echo writing codeml ctl file for $aln
	touch $ctldir/$abbrev.ctl
	echo "seqfile = $aln * sequence data filename" >> $ctldir/$abbrev.ctl
	echo "outfile = $codemldir/$iname.txt * main result file name" >> $ctldir/$abbrev.ctl
	echo "*treefile = cluster_1_tree2.nw" >> $ctldir/$abbrev.ctl
	echo "noisy = 0 * 0,1,2,3,9: how much rubbish on the screen" >> $ctldir/$abbrev.ctl
	echo "verbose = 0 * 1:detailed output" >> $ctldir/$abbrev.ctl
	echo "runmode = -2 * -2:pairwise" >> $ctldir/$abbrev.ctl
	echo "seqtype = 1 * 1:codons" >> $ctldir/$abbrev.ctl
	echo "CodonFreq = 2 * 0:equal, 1:F1X4, 2:F3X4, 3:F61" >> $ctldir/$abbrev.ctl
	echo "model = 1 *" >> $ctldir/$abbrev.ctl
	echo "NSsites = 0 *" >> $ctldir/$abbrev.ctl
	echo "icode = 0 * 0:universal code" >> $ctldir/$abbrev.ctl
	echo "fix_kappa = 1 * 1:kappa fixed, 0:kappa to be estimated" >> $ctldir/$abbrev.ctl
	echo "kappa = 1 * initial or fixed kappa" >> $ctldir/$abbrev.ctl
	echo "fix_omega = 0 * 1:omega fixed, 0:omega to be estimated" >> $ctldir/$abbrev.ctl
	echo "omega = 0.5 * initial omega value" >> $ctldir/$abbrev.ctl
	echo "*ndata = 1" >> $ctldir/$abbrev.ctl
done

