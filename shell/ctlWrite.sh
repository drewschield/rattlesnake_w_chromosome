dir=$1
ctldir=$2
codemldir=$3
for aln in $dir/*.pal2nal; do
	name=$aln
	iname=${name%.pal2nal*}		# removes pal2nal from name
	iname=${iname##*/}			# removes directory name
	echo writing codeml ctl file for $aln
	touch $ctldir/$iname.ctl
	echo "seqfile = $aln * sequence data filename" >> $ctldir/$iname.ctl
	echo "outfile = $codemldir/$iname.txt * main result file name" >> $ctldir/$iname.ctl
	echo "*treefile = cluster_1_tree2.nw" >> $ctldir/$iname.ctl
	echo "noisy = 0 * 0,1,2,3,9: how much rubbish on the screen" >> $ctldir/$iname.ctl
	echo "verbose = 0 * 1:detailed output" >> $ctldir/$iname.ctl
	echo "runmode = -2 * -2:pairwise" >> $ctldir/$iname.ctl
	echo "seqtype = 1 * 1:codons" >> $ctldir/$iname.ctl
	echo "CodonFreq = 2 * 0:equal, 1:F1X4, 2:F3X4, 3:F61" >> $ctldir/$iname.ctl
	echo "model = 1 *" >> $ctldir/$iname.ctl
	echo "NSsites = 0 *" >> $ctldir/$iname.ctl
	echo "icode = 0 * 0:universal code" >> $ctldir/$iname.ctl
	echo "fix_kappa = 1 * 1:kappa fixed, 0:kappa to be estimated" >> $ctldir/$iname.ctl
	echo "kappa = 1 * initial or fixed kappa" >> $ctldir/$iname.ctl
	echo "fix_omega = 0 * 1:omega fixed, 0:omega to be estimated" >> $ctldir/$iname.ctl
	echo "omega = 0.5 * initial omega value" >> $ctldir/$iname.ctl
	echo "*ndata = 1" >> $ctldir/$iname.ctl
done

