bed=$1
feature=$2
refname=$3
for sample in CV0011 CV0629 CV0646 CV0650; do
	mosdepth -t 4 --fast-mode -n -b $bed -f ../CroVir_genome_L77pg_16Aug2017.final_rename.fasta ./results/female.$sample.$refname.$feature ./bam/$sample
done
