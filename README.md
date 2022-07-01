# Rattlesnake W chromosome analysis

![Rattlesnake W chromosome analysis](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps used to assemble, annotate, and characterize the female-specific W chromosome of the prairie rattlesnake (_Crotalus viridis_), including comparative genomic analyses with other caenophidian snake and amniote species. This workflow is a companion to the methods described in Schield et al. (in review).

Lists and reference files can be found in the `resources` directory. Shell and Python scripts are in `shell` and `python` directories, respectively. R scripts are in the `R` directory. Note that you may need to adjust the organization of path/file locations to suit your environment. This workflow assumes that you return to the main working directory after each major section.

Feel free to email me at drew.schield[at]colorado.edu with any questions.

## Contents

* [Software and dependencies](#software-and-dependencies)
* [Female genome assembly](#female-genome-assembly)
* [Identification of W chromosome scaffolds](#identification-of-w-chromosome-scaffolds)
* [W chromosome annotation](#w-chromosome-annotation)
* [ZW gametolog divergence](#zw-gametolog-divergence)
* [Comparative Z chromosome mapping in caenophidan snakes](#comparative-z-chromosome-mapping-in-caenophidian-snakes)
* [GC and CpG content](#gc-and-cpg-content)
* [Repeat content](#repeat-content)
* [mdg4 element GC-richness](#mdg4-element-gc-richness)
* [Annotation of full-length LTR elements](#annotation-of-full-length-ltr-elements)
* [Refugium index analysis](#refugium-index-analysis)
* [Toxicity index analysis](#toxicity-index-analysis)
* [Gene expression analysis](#gene-expression-analysis)
* [Gene decay on the W chromosome](#gene-decay-on-the-w-chromosome)
* [W-specific gene duplications](#w-specific-gene-duplications)
* [Sex-linked divergence between pitvipers](#sex-linked-divergence-between-pitvipers)
* [ZW gametolog GC3 analysis](#zw-gametolog-GC3-analysis)
* [Appendix: Indian cobra analysis](#appendix-1-indian-cobra-analysis) 

## Software and dependencies

The analysis sections below use the following software and dependencies and assume they are on the user path:

* [FastQC](https://github.com/s-andrews/FastQC)
* [MultiQC](https://multiqc.info/)
* [10x Genomics Supernova](https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/running) (v2.1.1)
* [NCBI BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
* [MashMap](https://github.com/marbl/MashMap)
* [SRA toolkit](https://github.com/ncbi/sra-tools)
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
* [Maker](https://www.yandell-lab.org/software/maker.html)
* [Agouti](https://github.com/svm-zhang/AGOUTI)
* [Repeatmasker](https://www.repeatmasker.org/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bwa](http://bio-bwa.sourceforge.net/)
* [htslib](http://www.htslib.org/)
* [samtools](http://www.htslib.org/)
* [bgzip](http://www.htslib.org/)
* [tabix](http://www.htslib.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [mosdepth](https://github.com/brentp/mosdepth)
* [GffRead](https://github.com/gpertea/gffread)
* [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
* [STAR](https://github.com/alexdobin/STAR)
* [R](https://cran.r-project.org/)

Note, I installed a number of these programs to my [conda](https://docs.conda.io/en/latest/) environment.

## Female genome assembly

### 10x Genomics Chromium linked-read sequencing

We generated a 10x Genomics Chromium library for a female prairie rattlesnake from snap-frozen liver tissue. The library was sequenced on an Illumina NovaSeq 6000 using 150 bp paired-end reads.

### FastQC/MultiQC analysis

Assess read quality using FastQC and summarize using MultiQC.

#### Set up environment

```
mkdir fastq
mkdir fastqc
mkdir genome_crotalus_female
```

Raw linked-read data for the female prairie rattlesnake should be placed in the `fastq` directory.

#### Run FastQC analysis

```
cd fastqc
fastqc --threads 16 -o . ../fastq/PitViper_042019_S2_L004_R1_001.fastq.gz ../fastq/PitViper_042019_S2_L004_R2_001.fastq.gz
```

#### Summarize output with MultiQC

```
multiqc .
```

The output `multiqc_report.html` can be viewed and data/plots can be exported.

#### Examine results in R

Run `multiqc_summary.R` for a visual summary of FastQC results using summary files in the `resources` directory.

### Assembly

Use Supernova to assembly linked-read data into pseudohaplotype scaffolds.

The assembler is optimized to work with 56-fold coverage, requiring a number of input reads based on genome size and this coverage threshold. Using an estimated 1.5 Gb genome size for prairie rattlesnake, read lengths (150 bp) and 56-fold coverage, the number of input reads is as follows: <br />

1,500,000,000 x 56 / 150 = __560,000,000 input reads__.

```
supernova run --maxreads 560000000 --id=CV0650_female_Cviridis --fastqs=./fastq/
```

Generate pseudohaplotype fasta assembly with a minimum scaffold length of 10 kb.

```
supernova mkoutput --style=pseudohap2 --asmdir=./CV0650_female_Cviridis/outs/assembly/ --outprefix=CV0650_10xAssembly_Round3 --minsize=10000 --headers=short 
```

Move assembly to female genome directory and clean up intermediates.

```
mv ./CV0650_female_Cviridis/outs/assembly/CV0650_10xAssembly_Round3_pseudohap.*.fasta ./genome_crotalus_female
rm -r ./CV0650_female_Cviridis/
```

## Identification of W chromosome scaffolds

Use homology with the male prairie rattlesnake reference genome and relative female/male read depths on female scaffolds to identify candidate W chromosome scaffolds.

### Set up environment

```
mkdir W_chromosome_identification
mkdir genome_crotalus
```

Retrieve the [assembly](https://figshare.com/ndownloader/files/16522091), [gene annotation](https://figshare.com/ndownloader/files/16522322), and [repeat annotation](https://figshare.com/ndownloader/files/16522487) for the male prairie rattlesnake reference. The assembly is also available from [NCBI]().

```
wget https://figshare.com/ndownloader/files/16522091 ./genome_crotalus/CroVir_genome_L77pg_16Aug2017.final_rename.fasta.gz
wget https://figshare.com/ndownloader/files/16522322 ./genome_crotalus/CroVir_rnd1.all.maker.final.homologIDs.gff.gz 
wget https://figshare.com/ndownloader/files/16522487 ./genome_crotalus/CroVir_genome_L77pg_16Aug2017.repeat.masked.final.out.gff3.gz
```

### Extract candidate W chromosome scaffolds using homology

W-linked sequences should not have stringent high-similarity hits to autosomes. Instead, these should either have hits to the Z chromosome or not hit anything in the male reference.

Use mashmap to compare the female and male assemblies.

```
mashmap -r ./genome_crotalus/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ./genome_crotalus/female/CV0650_10xAssembly_Round3_pseudohap.1.fasta -t 6 --perc_identity 95 -f one-to-one -o ./W_chromosome_identification/10xPsdo1_toCvv_pi95_1to1.mashmap.out
mashmap -r ./genome_crotalus/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ./genome_crotalus/female/CV0650_10xAssembly_Round3_pseudohap.2.fasta -t 6 --perc_identity 95 -f one-to-one -o ./W_chromosome_identification/10xPsdo2_toCvv_pi95_1to1.mashmap.out
```

From this can be parsed a list of female scaffolds meeting these expectations for W-linkage:
* Not homologous to autosomes
* Not homologous to the pseudoautosomal region (i.e., scaffold-Z:106800000-113984505)
* Homology either to 1) the Z chromosome or 2) no hits to the male reference

These lists can be found in `resources/CvvPseudo1_NoAutoHits_scaffIDs_02.27.20.txt` and `resources/CvvPseudo2_NoAutoHits_scaffIDs_02.27.20.txt`.

### Use relative read depths for female and male to confirm W-linkage

With a list of candidates based on homology, compare normalized read depths for a female and male across these scaffolds to find those with ratios expected for the female-specific W chromosome.

#### 1. Compare alternative methods for identifying sex-linked sequences

Different statistical thresholds for deciding whether a scaffold is W-linked or not may have different power and false positive rates. Two common approaches are comparing the log2 ratio of normalized female:male read depth (log2FM) and comparing the proportion of reads mapping to a sequence from males versus females. Here, we can use known W-linked sequence from the chicken (_Gallus gallus_) to ground truth these thresholds and compare.

Thresholds to compare are:
* log2FM > 1
* female mapping proportion > Q3 + 1.5*IQR

Retrieve chicken reference from NCBI and read data from SRA (accessions SRR958465 [male] and SRR958466 [female]).

```
mkdir genome_gallus
cd genome_gallus
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCA_000002315.5_GRCg6a/GCA_000002315.5_GRCg6a_genomic.fna.gz
mkdir fastq
cd fastq
fastq-dump --split-files --gzip SRR958465
fastq-dump --split-files --gzip SRR958466
mv SRR958465_1.fastq.gz gallus_male_SRR958465_1.fastq.gz
mv SRR958465_2.fastq.gz gallus_male_SRR958465_2.fastq.gz
mv SRR958466_1.fastq.gz gallus_female_SRR958466_1.fastq.gz
mv SRR958466_2.fastq.gz gallus_female_SRR958466_2.fastq.gz
cd ..
```

Index reference, map reads, and index outputs.

```
bwa index ./genome_gallus/GCA_000002315.5_GRCg6a_genomic.fna
mkdir ./W_chromosome_identification/comparative_W_coverage/gallus
mkdir ./W_chromosome_identification/comparative_W_coverage/gallus/bam
bwa mem -t 16 ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.fna ./genome_gallus/fastq/gallus_male_SRR958465_1.fastq.gz ./genome_gallus/fastq/gallus_male_SRR958465_2.fastq.gz | samtools sort -O bam -T harp -o ./W_chromosome_identification/comparative_W_coverage/gallus/bam/gallus_male_GRCg6a.bam -
bwa mem -t 16 ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.fna ./genome_gallus/fastq/gallus_female_SRR958466_1.fastq.gz ./genome_gallus/fastq/gallus_female_SRR958466_2.fastq.gz | samtools sort -O bam -T inni -o ./W_chromosome_identification/comparative_W_coverage/gallus/bam/gallus_female_GRCg6a.bam -
samtools index ./W_chromosome_identification/comparative_W_coverage/gallus/bam/gallus_female_GRCg6a.bam 
samtools index ./W_chromosome_identification/comparative_W_coverage/gallus/bam/gallus_male_GRCg6a.bam 
```

Determine reasonable window size for analysis in chicken.

The chicken scaffolds are highly contiguous. In order to make a reasonable comparison to the prairie rattlesnake female scaffolds, we can determine the mean length of female rattlesnake scaffolds, and make a bed file for the chicken genome with window sizes equal to this length.

The mean female rattlesnake scaffold length is __10013.8 bp__, so setting chicken analysis to 10 kb windows seems reasonable.

Make a bed file of chicken scaffolds in 10 kb windows after extracting 'genome' file using the script `fastq_seq_length.py`.

```
python fasta_seq_length.py ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.fna > ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.scaffold_lengths.txt
bedtools makewindows -g ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.scaffold_lengths.txt -w 10000 > ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.10kb.window.bed
```

Calculate mean depth per 10 kb window using mosdepth.

```
mosdepth -t 4 --fast-mode -n -b ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.10kb.window.bed ./W_chromosome_identification/comparative_W_coverage/gallus/gallus_female_GRCg6a.10kb ./W_chromosome_identification/comparative_W_coverage/gallus/bam/gallus_female_GRCg6a.bam
mosdepth -t 4 --fast-mode -n -b ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.10kb.window.bed ./W_chromosome_identification/comparative_W_coverage/gallus/gallus_male_GRCg6a.10kb ./W_chromosome_identification/comparative_W_coverage/gallus/bam/gallus_male_GRCg6a.bam
```

Examine results in R.

Run `W_identification_power.R` to calculate power to detect W-linked sequence (and false positive rate) in the chicken using the log2FM and Q3 + 1.5*IQR methods.

These analyses demonstrate that the __log2FM >= 1__ threshold method has high power and low false positive rate, and so will be used in the identification of candidate W scaffolds in prairie rattlesnake.

#### 2. Calculate relative female:male coverage in prairie rattlesnake

Set up environment.

```
mkdir W_chromosome_identification/coverage_exp_crotalus
mkdir W_chromosome_identification/coverage_exp_crotalus/fastq
mkdir W_chromosome_identification/coverage_exp_crotalus/bam
```

Retrieve read data for male and female prairie rattlesnakes.

```
./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R1_001.fastq.gz
./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R1_001.fastq.gz
./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R1_001.fastq.gz
./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R1_001.fastq.gz
```

Index pseudohaplotype reference fasta files.

```
bwa index ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.1.fasta
bwa index ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.2.fasta
```

Map to female reference.

```
bwa mem -t 16 ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.1.fasta ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R1_001.fastq.gz ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R2_001.fastq.gz | samtools sort -O bam -T bimb -o ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0007_CV0650_pseudohap1.bam -
bwa mem -t 16 ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.2.fasta ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R1_001.fastq.gz ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R2_001.fastq.gz | samtools sort -O bam -T chup -o ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0007_CV0650_pseudohap2.bam -
bwa mem -t 16 ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.1.fasta ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R1_001.fastq.gz ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R2_001.fastq.gz | samtools sort -O bam -T flug -o ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0011_CV0650_pseudohap1.bam -
bwa mem -t 16 ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.2.fasta ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R1_001.fastq.gz ./W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R2_001.fastq.gz | samtools sort -O bam -T posh -o ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0011_CV0650_pseudohap2.bam -
```

Index output bam files.

```
samtools index ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0011_CV0650_pseudohap1.bam
samtools index ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0011_CV0650_pseudohap2.bam
samtools index ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0007_CV0650_pseudohap1.bam
samtools index ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0007_CV0650_pseudohap2.bam
```

Calculate mean depth per scaffold using mosdepth.

```
mosdepth -t 4 --fast-mode -n ./W_chromosome_identification/coverage_exp_crotalus/CV0011_CV0650_pseudohap1 ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0011_CV0650_pseudohap1.bam
mosdepth -t 4 --fast-mode -n ./W_chromosome_identification/coverage_exp_crotalus/CV0011_CV0650_pseudohap2 ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0011_CV0650_pseudohap2.bam
mosdepth -t 4 --fast-mode -n ./W_chromosome_identification/coverage_exp_crotalus/CV0007_CV0650_pseudohap1 ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0007_CV0650_pseudohap1.bam
mosdepth -t 4 --fast-mode -n ./W_chromosome_identification/coverage_exp_crotalus/CV0007_CV0650_pseudohap2 ./W_chromosome_identification/coverage_exp_crotalus/bam/CV0007_CV0650_pseudohap2.bam
```

Examine results and identify candidate W chromosome scaffolds in R.

Run `W_identification_crotalus.R` to compare normalized female:male read depths, parse scaffolds with log2FM > 1, and cross-reference with candidate W scaffolds from homology search.

#### 3. Filter and parse candidate W chromosome scaffolds

Filter scaffolds with no hits to autosomes.

```
cd ./W_chromosome_identification/
touch pseudohap1_no-auto_log2FM_thresh_scaffolds.txt; for i in `cat CvvPseudo1_NoAutoHits_scaffIDs_02.27.20.txt`; do awk 'BEGIN{OFS='\t'}{if ($1=='$i') print $0}' pseudohap1_log2FM_thresh_scaffolds.txt >> pseudohap1_no-auto_log2FM_thresh_scaffolds.txt; done
touch pseudohap2_no-auto_log2FM_thresh_scaffolds.txt; for i in `cat CvvPseudo2_NoAutoHits_scaffIDs_02.27.20.txt`; do awk 'BEGIN{OFS='\t'}{if ($1=='$i') print $0}' pseudohap2_log2FM_thresh_scaffolds.txt >> pseudohap2_no-auto_log2FM_thresh_scaffolds.txt; done
cd ..
```

The final lists of candidate W scaffolds are in `CANDIDATE_W_pseudohap1_scaffold.list.txt` and `CANDIDATE_W_pseudohap1_scaffold.list.txt` in the `resources` directory.

Extract sequence for candidate W scaffolds using `scaffold_list_extractor.py`.

```
mkdir ./W_chromosome_identification/candidate_W/
python scaffold_list_extractor.py ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.1.fasta ./W_chromosome_identification/CANDIDATE_W_pseudohap1_scaffold.list.txt ./W_chromosome_identification/candidate_W/pseudohaplotype1.candidate_W.fasta
python scaffold_list_extractor.py ./genome_crotalus_female/CV0650_10xAssembly_Round3_pseudohap.2.fasta ./W_chromosome_identification/CANDIDATE_W_pseudohap2_scaffold.list.txt ./W_chromosome_identification/candidate_W/pseudohaplotype2.candidate_W.fasta
cat ./W_chromosome_identification/candidate_W/pseudohaplotype1.candidate_W.fasta ./W_chromosome_identification/candidate_W/pseudohaplotype2.candidate_W.fasta > ./W_chromosome_identification/candidate_W/Cviridis_CV0650_candidate_W.fasta
```

## W chromosome annotation

### Repeat elements

Genomic repeats on the W chromosome were annotated using RepeatMasker leveraging the following databases:
* Bov-B/CR1 LINE ([Pasquesi et al. 2018](https://www.nature.com/articles/s41467-018-05279-1))
* RepBase release 20181026 ([Bao et al. 2015](https://link.springer.com/article/10.1186/s13100-015-0041-9))
* Snake-specific library ([Schield et al. 2019](https://genome.cshlp.org/content/29/4/590))

### Gene annotation and transcript-based scaffold improvement

Annotation of protein-coding genes included empirical evidence from a transcriptome assembly for prairie rattlesnake from [Schield et al. (2019)](https://genome.cshlp.org/content/29/4/590) and a novel transcriptome based on RNAseq data from female prairie rattlesnake tissues generated using Trinity. Additional empirical evidence came from protein datasets for [*Anolis carolinensis*](https://www.nature.com/articles/nature10390), [*Python molurus bivittatus*](https://www.pnas.org/doi/abs/10.1073/pnas.1314475110), [*Thamnophis sirtalis*](https://academic.oup.com/gbe/article/10/8/2110/5061318?login=true), [*Ophiophagus hannah*](https://www.pnas.org/doi/abs/10.1073/pnas.1314702110), and [*Deinagkistrodon acutus*](https://www.nature.com/articles/ncomms13107).

Empirical evidence was used in an initial annotation using MAKER using default settings, except for specifying `max_dna_len=300000` and `split_hit=20000`. The results were used to optimize gene prediction parameters in Augustus, then passed to a second run of MAKER.

The *de novo* transcriptome and candidate W assembly data were used to improve contiguity of W chromosome scaffolds using Agouti.

Repeat and gene annotations were repeated following scaffold improvement.

The W chromosome assembly and annotation files are in the `resources/annotation/` directory. The assembly is also deposited at [NCBI](xxxx).

## ZW gametolog divergence

We can investigate divergence between ZW gametologs and the presence of evolutionary strata based on spatial clustering of gametolog pairs with similar divergence along the Z chromosome.

### 1. Estimation of a lineage-specific mutation rate

First, use divergence between rattlesnake and anole lizard 1:1 orthologs and their known divergence time to calculate a mutation rate estimate to scale divergence estimates between ZW gametologs.

#### Set up environment

```
mkdir divergence_crotalus_anolis
mkdir divergence_crotalus_anolis/ortholog_sequences
mkdir divergence_crotalus_anolis/ortholog_alignments
mkdir divergence_crotalus_anolis/ctl
mkdir divergence_crotalus_anolis/codeml
```

#### Retrieve anole data

```
mkdir genome_anolis
cd genome_anolis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/745/GCF_000090745.1_AnoCar2.0/GCF_000090745.1_AnoCar2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/745/GCF_000090745.1_AnoCar2.0/GCF_000090745.1_AnoCar2.0_genomic.gff.gz
cd ..
```

#### Extract CDS sequences

```
gffread -x ./divergence_crotalus_anolis/anolis.cds.fasta -g ./genome_anolis/GCF_000090745.1_AnoCar2.0_genomic.fna ./genome_anolis/GCF_000090745.1_AnoCar2.0_genomic.gff
grep -v '#' ./genome_crotalus/CroVir_rnd1.all.maker.final.homologIDs.gff | gffread -x ./divergence_crotalus_anolis/crotalus.cds.fasta -g ./genome_crotalus/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -
```

The initial grep command for *Crotalus* removes the commented entries in the GFF (GffRead doesn't want to see these).

#### Identify 1:1 orthologs using tBLASTx

Make BLAST databases for reciprocal searches.

```
makeblastdb -dbtype nucl -in ./divergence_crotalus_anolis/anolis.cds.fasta
makeblastdb -dbtype nucl -in ./divergence_crotalus_anolis/crotalus.cds.fasta
```

Perform reciprocal tBLASTx searches.

```
tblastx -num_threads 8 -max_target_seqs 5 -max_hsps 1 -evalue 0.00001 -outfmt "6 qacc sacc evalue bitscore qstart qend sstart send" -db ./divergence_crotalus_anolis/crotalus.cds.fasta -query ./divergence_crotalus_anolis/anolis.cds.fasta -out ./divergence_crotalus_anolis/tblastx_anolis2crotalus.cds.txt
tblastx -num_threads 4 -max_target_seqs 5 -max_hsps 1 -evalue 0.00001 -outfmt "6 qacc sacc evalue bitscore qstart qend sstart send" -db ./divergence_crotalus_anolis/anolis.cds.fasta -query ./divergence_crotalus_anolis/crotalus.cds.fasta -out ./divergence_crotalus_anolis/tblastx_crotalus2anolis.cds.txt
```

Extract reciprocal best BLAST hits using `RBH_comma.py` script, adapted from the script written by [Daren Card](https://github.com/darencard). 

```
python RBH_comma.py ./divergence_crotalus_anolis/tblastx_anolis2crotalus.cds.txt ./divergence_crotalus_anolis/tblastx_crotalus2anolis.cds.txt ./divergence_crotalus_anolis/orthologs_crotalus_anolis.one2one.txt
```

This identified 12,367 ortholog pairs.

Remove pairs that are Z-linked in rattlesnake.

```
grep -v 'scaffold-Z' ./divergence_crotalus_anolis/orthologs_crotalus_anolis.one2one.txt > ./divergence_crotalus_anolis/orthologs_crotalus_anolis.one2one.autosome.txt
```

This keeps __11,277 autosomal ortholog pairs__.

Extract fasta sequences for orthologs using `makeOrthoFasta.py`.

```
cd divergence_crotalus_anolis
python makeOrthoFasta.py orthologs_crotalus_anolis.one2one.autosome.txt crotalus.cds.fasta anolis.cds.fasta
```

This writes sequence pairs to `./divergence_crotalus_anolis/ortholog_sequences/`.

Run `translate_fasta.py` to translate nucleotide fastas to amino acid sequences.

```
for fasta in ./ortholog_sequences/*.fna; do python ./python/translate_fasta.py $fasta; done
```

#### Align 1:1 orthologs using Clustal Omega

Run `alignClustal.sh` to align amino acid sequences.

```
sh alignClustal.sh ./ortholog_sequences ./ortholog_alignments
```

Generate codon-aware nucleotide alignments using PAL2NAL with `convertPAL2NAL.sh` script.

```
sh convertPAL2NAL.sh ./ortholog_alignments ./ortholog_sequences
```

Write a control file for PAML for each alignment using `ctlWrite.sh`.

```
sh ctlWrite.sh ./ortholog_alignments ./ctl ./codeml
```

PAML barfs on the long-format detailed names. These steps will ensure each of the files has simple 'crotalus' and 'anolis' sequence headers. The second sed command takes care of the anole entries because they all have '-' in them and the rattlesnake lines have already been overwritten.

```
for fasta in ./ortholog_sequences/*.fna; do sed -i '/-scaffold-/c\>crotalus' $fasta; sed -i '/-/c\>anolis' $fasta; done
for fasta in ./ortholog_sequences/*.faa; do sed -i '/-scaffold-/c\>crotalus' $fasta; sed -i '/-/c\>anolis' $fasta; done
for fasta in ./ortholog_alignments/*.pal2nal; do sed -i '/-scaffold-/c\crotalus' $fasta; sed -i '/-/c\anolis' $fasta; done
```

The aligned amino acids need a bit more attention, since gaps are denoted as '-':

```
for fasta in ./ortholog_alignments/*.faa; do sed -i '/-scaffold-/c\>crotalus' $fasta; sed -i '/_/c\>anolis' $fasta; done
for fasta in ./ortholog_alignments/*.faa; do sed -i '/id-LOC/c\>anolis' $fasta; done
```

Check that headers have all been reformatted:

```
for i in ./ortholog_sequences/*.fna; do grep '>anolis' $i; done | wc -l
for i in ./ortholog_sequences/*.faa; do grep '>anolis' $i; done | wc -l
for i in ./ortholog_alignments/*.faa; do grep '>anolis' $i; done | wc -l	
for i in ./ortholog_alignments/*.pal2nal; do grep 'anolis' $i; done | wc -l
```

All should equal 11,277.

#### Calculate divergence statistics using CODEML

Run CODEML on each alignment, calling the respective control file.

```
for control in ./ctl/*.ctl; do codeml $control; done
```

Run `parseCodeml.sh` to write a table of divergence statistics for all of the alignments. This calls the `parse_codeml_output.py` script available [here](https://github.com/faylward/dnds), modified to print 'NA's for filtered results.

```
sh parseCodeml.sh > crotalus_anolis_ortholog.autosome.dnds.txt
```

#### Convert of rattlesnake-anole divergence estimates to time

Run `divergence_crotalus_anolis.R` to convert divergence estimates to time, including tranformation to account sex-linked mutation rates.

An alternative generalized squamate mutation rate from four-fold degenerate sites is available from [Green et al. (2014)](https://www.science.org/doi/full/10.1126/science.1254449).

### 2. Identification of 1:1 ZW gametologs in prairie rattlesnake

Identify 1:1 gametologs on the Z and W chromosomes to estimate divergence.

#### Set up environment

```
mkdir divergence
mkdir divergence/crotalus/
cd ./divergence/crotalus/
```

#### Extract Z and W chromosome CDS sequences

```
grep 'scaffold-Z' ../../genome_crotalus/CroVir_rnd1.all.maker.final.homologIDs.gff | gffread -x chrZ.cds.fasta -g ../../genome_crotalus/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -
gffread -x chrW.cds.fasta -g ../../resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta ../../resources/annotation/croVir_Wscaff_rnd2_alt.all.maker.noseq.gff3
```

#### Make BLAST databases for CDS sequences

```
makeblastdb -dbtype nucl -in chrZ.cds.fasta
makeblastdb -dbtype nucl -in chrW.cds.fasta
```

#### Perform reciprocal tBLASTx searches

```
tblastx -num_threads 8 -max_hsps 1 -evalue 0.00001 -outfmt "6 qacc sacc evalue bitscore qstart qend sstart send" -db chrW.cds.fasta -query chrZ.cds.fasta -out tblastx_Z2W.cds.txt
tblastx -num_threads 8 -max_hsps 1 -evalue 0.00001 -outfmt "6 qacc sacc evalue bitscore qstart qend sstart send" -db chrZ.cds.fasta -query chrW.cds.fasta -out tblastx_W2Z.cds.txt
```

#### Process results using `RBH_comma.py` script

```
python RBH_comma.py tblastx_Z2W.cds.txt tblastx_W2Z.cds.txt gametologs.cds.one2one.txt
```

Format a tab-delimited version to paste together with gene coordinate and annotation details.

```
sed 's/,/\t/g' gametologs.cds.one2one.txt > gametologs.cds.one2one.fix.txt
```

#### Obtain Z chromosome locations for 1:1 gametologs

Coordinates are needed for ZW gametolog pairs to interpret divergence estimates in the context of the Z chromosome.

Extract gene IDs from the Z chromosome CDS.

```
grep '>' chrZ.cds.fasta | cut -d'>' -f2 > cds.Z_geneID.list
```

Query GFF for coordinates of each gene.

```
touch cds.Z_geneID_coordinates.txt; for mrna in `cat cds.Z_geneID.list`; do grep -w "$mrna" ../../genome_crotalus/CroVir_rnd1.all.maker.final.homologIDs.gff | awk 'BEGIN{OFS="\t"} $3 == "mRNA" {print $1,$4,$5,$9}' - >> cds.Z_geneID_coordinates.txt; done
```

Note that the fourth column also contain the transcript IDs.

Extract Z chromosome coordinates for ZW gametolog pairs.

```
touch ZW_homology.Z_geneID_coordinates.txt; cat gametologs.cds.one2one.txt | while read line; do chrz=`echo $line | cut -d, -f2`; grep -w $chrz cds.Z_geneID_coordinates.txt >> ZW_homology.Z_geneID_coordinates.txt; done
```

Paste coordinates for ZW gametologs to annotation details.

```
pr -mts$'\t' <(cut -f1,2,3 ZW_homology.Z_geneID_coordinates.txt) gametologs.cds.one2one.fix.txt <(cut -f4 ZW_homology.Z_geneID_coordinates.txt) > gametologs.cds.one2one.Z_geneID_coordinates.txt
```

The columns of the output are in the following order:

chromosome, start position, end position, W transcript, Z transcript, annotation details

### 3. Alignment of ZW gametologs

#### Set up environment

```
cd ./divergence/crotalus/
mkdir ZW_gametolog_seq
mkdir ZW_gametolog_aln
```

#### Extract ZW gametolog sequence files

Run `make_ZW_gametolog_fasta.py` to write fasta file per ZW gametolog pair.

```
python make_ZW_gametolog_fasta.py gametologs.cds.one2one.Z_geneID_coordinates.txt chrZ.cds.fasta chrW.cds.fasta
```

This writes the output files to the `ZW_gametolog_seq` directory.

#### Translate sequences to amino acid

Run `translate_gametolog_fasta.py` on input nucleotide fasta files.

```
for fasta in ZW_gametolog_seq/*.fna; do python translate_gametolog_fasta.py $fasta; done
```

#### Align ZW gametologs

First, run `alignClustal_gametologs.sh` to align amino acid sequences using Clustal Omega.

```
sh alignClustal_gametologs.sh ZW_gametolog_seq ZW_gametolog_aln
```

Then run `convertPAL2NAL_gametologs.sh` to convert to codon-based nucleotide alignments.

```
sh convertPAL2NAL_gametologs.sh ZW_gametolog_aln ZW_gametolog_seq
```

### 4. Analysis of ZW gametologs in CODEML

#### Set up environment

```
mkdir ctl
mkdir codeml
```

#### Generate CODEML control files

Run `ctlWrite_gametologs.sh` to generate control files for analysis of each alignment, specifying the directories with alignments, control files, and where the results of codeml will be written.

```
sh ctlWrite_gametologs.sh ZW_gametolog_aln ctl codeml
```

#### Perform CODEML analyses

Run CODEML to calculate divergence statistics per alignment.

```
for control in ./ctl/*.ctl; do codeml $control; done
```

Parse CODEML results using `parseCodeml_gametologs.sh`.

```
sh parseCodeml_gametologs.sh > crotalusZW_gametolog.dnds.txt
```

## Comparative Z chromosome mapping in caenophidian snakes

## GC and CpG content

## Repeat content

## mdg4 retroelement GC-richness

## Annotation of full-length LTR elements

## Refugium index analysis

## Toxicity index analysis

## Gene expression analysis

## Gene decay on the W chromosome

## W-specific gene duplications

## Sex-linked divergence between pitvipers

## ZW gametolog GC3 analysis

## Appendix: Indian cobra analysis










































