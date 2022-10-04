# Rattlesnake W chromosome analysis

![Rattlesnake W chromosome analysis](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps used to assemble, annotate, and characterize the female-specific W chromosome of the prairie rattlesnake (_Crotalus viridis_), including comparative genomic analyses with other caenophidian snake and amniote species. This workflow is a companion to the methods described in Schield et al. (_in review_).

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
* [mdg4 retroelement GC-richness](#mdg4-retroelement-gc-richness)
* [Annotation of full-length LTR elements](#annotation-of-full-length-ltr-elements)
* [Refugium & toxicity index analyses](#refugium-toxicity-index-analyses)
* [Gene expression analysis](#gene-expression-analysis)
* [W-specific gene duplications](#w-specific-gene-duplications)
* [Sex-linked divergence between pitvipers](#sex-linked-divergence-between-pitvipers)
* [ZW gametolog GC3 analysis](#zw-gametolog-GC3-analysis)
* [Appendix: Indian cobra analysis](#appendix-indian-cobra-analysis) 

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
* [GenomeTools](https://github.com/genometools/genometools)
* [HMMER](https://github.com/EddyRivasLab/hmmer)
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

### De novo transcriptome assembly
We assembled a *de novo* transcriptome using RNA-seq data from 18 female tissue samples (See Supplementary Table S2).

#### Subsampling reads with seqtk
Reads were randomly subsampled using [seqtk](www.github.com/lh3/seqtk) to retain 2 million reads per sample. 
```bash
# Example command:
seqtk sample -s100 Hrt_24h_CV8_Cvv11_GGCTAC__R1_merged.fastq.gz 2000000 > subsample/Hrt_24h_CV8_Cvv11_GGCTAC__R1_merged.subsamp.fastq
```

#### Transcriptome assembly with Trinity
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq) was used to assemble a transcriptome using the --trimmomatic flag to incorporate upfront quality trimming and otherwise default parameters.
```bash
./trinityrnaseq_r20140717/Trinity --seqType fq --JM 70G --trimmomatic --CPU 12 --SS_lib_type FR --left ./raw_reads/subsample/Hrt_24h_CV8_Cvv11_GGCTAC__R1_merged.subsamp.fastq ./raw_reads/subsample/Hrt_96h_CV10_Cvv15_CCGTCC__R1_merged.subsamp.fastq ./raw_reads/subsample/Hrt_96h_CV11_Cvv16_GTAGAG__R1_merged.subsamp.fastq ./raw_reads/subsample/Hrt_96h_CV9_Cvv14_ATGTCA__R1_merged.subsamp.fastq ./raw_reads/subsample/Kid_24h_CV3_Cvv19_CAGGCG_L005_R1_001.subsamp.fastq ./raw_reads/subsample/Kid_24h_CV8_Cvv17_CACGAT_L005_R1_001.subsamp.fastq ./raw_reads/subsample/Kid_96h_CV10_Cvv24_CTAGCT_L005_R1_001.subsamp.fastq ./raw_reads/subsample/Kid_96h_CV11_Cvv21_CATTTT_L005_R1_001.subsamp.fastq ./raw_reads/subsample/Kid_96h_CV9_Cvv23_CGGAAT_L005_R1_001.subsamp.fastq ./raw_reads/subsample/liver_24h_CV3_Cvv08_CGTACG_L005_R1_001.subsamp.fastq ./raw_reads/subsample/Liver_24h_CV8_Cvv05_GTGAAA_L005_R1_001.subsamp.fastq ./raw_reads/subsample/liver_96h_CV10_Cvv11_ACTGAT_L005_R1_001.subsamp.fastq ./raw_reads/subsample/liver_96h_CV11_Cvv10_GGTAGC_L005_R1_001.subsamp.fastq ./raw_reads/subsample/Liver_96h_CV9_Cvv12_ATGAGC_L005_R1_001.subsamp.fastq ./raw_reads/subsample/SI_24h_CV8_Cvv05_ACAGTG__R1_merged.subsamp.fastq ./raw_reads/subsample/SI_96h_CV11_Cvv08_ACTTGA__R1_merged.subsamp.fastq ./raw_reads/subsample/SI_96h_CV9_Cvv07_CAGATC__R1_merged.subsamp.fastq --right ./raw_reads/subsample/Hrt_24h_CV8_Cvv11_GGCTAC__R2_merged.subsamp.fastq ./raw_reads/subsample/Hrt_96h_CV10_Cvv15_CCGTCC__R2_merged.subsamp.fastq ./raw_reads/subsample/Hrt_96h_CV11_Cvv16_GTAGAG__R2_merged.subsamp.fastq ./raw_reads/subsample/Hrt_96h_CV9_Cvv14_ATGTCA__R2_merged.subsamp.fastq ./raw_reads/subsample/Kid_24h_CV3_Cvv19_CAGGCG_L005_R2_001.subsamp.fastq ./raw_reads/subsample/Kid_24h_CV8_Cvv17_CACGAT_L005_R2_001.subsamp.fastq ./raw_reads/subsample/Kid_96h_CV10_Cvv24_CTAGCT_L005_R2_001.subsamp.fastq ./raw_reads/subsample/Kid_96h_CV11_Cvv21_CATTTT_L005_R2_001.subsamp.fastq ./raw_reads/subsample/Kid_96h_CV9_Cvv23_CGGAAT_L005_R2_001.subsamp.fastq ./raw_reads/subsample/liver_24h_CV3_Cvv08_CGTACG_L005_R2_001.subsamp.fastq ./raw_reads/subsample/liver_24h_CV8_Cvv05_GTGAAA_L005_R2_001.subsamp.fastq ./raw_reads/subsample/liver_96h_CV10_Cvv11_ACTGAT_L005_R2_001.subsamp.fastq ./raw_reads/subsample/liver_96h_CV11_Cvv10_GGTAGC_L005_R2_001.subsamp.fastq ./raw_reads/subsample/liver_96h_CV9_Cvv12_ATGAGC_L005_R2_001.subsamp.fastq ./raw_reads/subsample/SI_24h_CV8_Cvv05_ACAGTG__R2_merged.subsamp.fastq ./raw_reads/subsample/SI_96h_CV11_Cvv08_ACTTGA__R2_merged.subsamp.fastq ./raw_reads/subsample/SI_96h_CV9_Cvv07_CAGATC__R2_merged.subsamp.fastq
```

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

Comparison of female and male read depths across the Z chromosome can reveal differences in relative timing of recombination suppression between the sex chromosomes if compared between different ZW snake species. This section contains details of the analysis comparing log2FM ratios in pitvipers to garter snake (_Thamnophis_).

### 1. Retrieve read data for snake species

These analyses will use read data for prairie rattlesnake described [above](#identification-of-w-chromosome-scaffolds), five-pace viper (_Deinagkistrodon acutus_) from [Yin et al. (2016)](https://www.nature.com/articles/ncomms13107), and pygmy rattlesnake (_Sistrurus miliarius_) and western terrestiral garter snake (_Thamnophis elegans_) from [Vicoso et al. (2013)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001643). 

#### Set up environment
```
mkdir coverage
mkdir coverage/fastq/
mkdir coverage/bam/
mkdir coverage/mosdepth_results/
cd coverage
```

#### Download read data for other species from NCBI SRA
```
fastq-dump --split-files --gzip ./fastq/SRR941260
fastq-dump --split-files --gzip ./fastq/SRR941259
fastq-dump --split-files --gzip ./fastq/SRR3223535
fastq-dump --split-files --gzip ./fastq/SRR3223527
fastq-dump --split-files --gzip ./fastq/SRR941263
fastq-dump --split-files --gzip ./fastq/SRR941261
```

Rename fastq files with species/sex identifiers.
```
mv ./fastq/SRR941260_1.fastq.gz ./fastq/thamnophis_mDNA_1.fastq.gz
mv ./fastq/SRR941260_2.fastq.gz ./fastq/thamnophis_mDNA_2.fastq.gz
mv ./fastq/SRR941259_1.fastq.gz ./fastq/thamnophis_fDNA_1.fastq.gz
mv ./fastq/SRR941259_2.fastq.gz ./fastq/thamnophis_fDNA_2.fastq.gz
mv ./fastq/SRR3223535_1.fastq.gz ./fastq/deinagkistrodon_mDNA_1.fastq.gz
mv ./fastq/SRR3223535_2.fastq.gz ./fastq/deinagkistrodon_mDNA_2.fastq.gz
mv ./fastq/SRR3223527_1.fastq.gz ./fastq/deinagkistrodon_fDNA_1.fastq.gz
mv ./fastq/SRR3223527_2.fastq.gz ./fastq/deinagkistrodon_fDNA_2.fastq.gz
mv ./fastq/SRR941263_1.fastq.gz ./fastq/sistrurus_mDNA_1.fastq.gz
mv ./fastq/SRR941263_2.fastq.gz ./fastq/sistrurus_mDNA_2.fastq.gz
mv ./fastq/SRR941261_1.fastq.gz ./fastq/sistrurus_fDNA_1.fastq.gz
mv ./fastq/SRR941261_2.fastq.gz ./fastq/sistrurus_fDNA_2.fastq.gz
```

#### Map reads to the rattlesnake genome
```
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./fastq/thamnophis_mDNA_1.fastq.gz ./fastq/thamnophis_mDNA_2.fastq.gz | samtools sort -O bam -T male -o ./bam/thamnophis_mDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./fastq/thamnophis_fDNA_1.fastq.gz ./fastq/thamnophis_fDNA_2.fastq.gz | samtools sort -O bam -T female -o ./bam/thamnophis_fDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./fastq/deinagkistrodon_mDNA_1.fastq.gz ./fastq/deinagkistrodon_mDNA_2.fastq.gz | samtools sort -O bam -T male -o ./bam/deinagkistrodon_mDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./fastq/deinagkistrodon_fDNA_1.fastq.gz ./fastq/deinagkistrodon_fDNA_2.fastq.gz | samtools sort -O bam -T female -o ./bam/deinagkistrodon_fDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./fastq/sistrurus_mDNA_1.fastq.gz ./fastq/sistrurus_mDNA_2.fastq.gz | samtools sort -O bam -T male -o ./bam/sistrurus_mDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./fastq/sistrurus_fDNA_1.fastq.gz ./fastq/sistrurus_fDNA_2.fastq.gz | samtools sort -O bam -T female -o ./bam/sistrurus_fDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ../W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R1_001.fastq.gz ../W_chromosome_identification/coverage_exp_crotalus/fastq/CV0007_S82_L004_R2_001.fastq.gz | samtools sort -O bam -T bimb -o ./bam/crotalus_mDNA.bam -
bwa mem -t 4 ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ../W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R1_001.fastq.gz ../W_chromosome_identification/coverage_exp_crotalus/fastq/CV0011_S70_L004_R2_001.fastq.gz | samtools sort -O bam -T flug -o ./bam/crotalus_fDNA.bam -
```

#### Filter mapping quality

Filter using samtools.
```
for map in ./bam/*.bam; do name=`echo $map | cut -d'/' -f3 | cut -d'.' -f1`; samtools view -q 30 -b $map > ./bam/$name.q30.bam; done
```

Index results.
```
for i in ./bam/*.q30.bam; do samtools index $i; done
```

#### Run mosdepth analysis

Perform analyses in 10 kb sliding windows using coordinates in `resources/CroVir_Dovetail_10kb_window.ChromAssigned.bed`.
```
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/thamnophis_mDNA ./bam/thamnophis_mDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/thamnophis_fDNA ./bam/thamnophis_fDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/deinagkistrodon_mDNA ./bam/deinagkistrodon_mDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/deinagkistrodon_fDNA ./bam/deinagkistrodon_fDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/sistrurus_mDNA ./bam/sistrurus_mDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/sistrurus_fDNA ./bam/sistrurus_fDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/crotalus_mDNA ./bam/crotalus_mDNA.bam
mosdepth -t 4 --fast-mode -n -b CroVir_Dovetail_10kb_window.ChromAssigned.bed -f ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta ./mosdepth_results/crotalus_fDNA ./bam/crotalus_fDNA.bam
```

Run `comparative_coverage_4species.R` to examine log2FM ratios across the Z chromosome in each species.

## GC and CpG content

With identified W-linked scaffolds, comparisons of GC and CpG content can be made between autosomal, Z-linked, and W-linked regions of the genome.

### 1. GC and CpG content in prairie rattlesnake

#### Set up environment
```
mkdir gc
cd gc
```

#### Calculate GC and CpG content per W chromosome scaffold

Run `scaffold_gc.py` to quantify GC content for each scaffold.
```
$python scaffold_gc.py ../resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta Cviridis_CV0650_candidate_W.rescaffold.rename.GC.txt
```

Run `scaffold_cpg.py` to quantify CpG content for each W chromosome scaffold.
```
python scaffold_cpg.py ../resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta Cviridis_CV0650_candidate_W.rescaffold.rename.CpG.txt
```

#### Calculate GC and CpG content in sliding windows on autosomes and Z chromosome

Run `slidingwindow_gc.py` to quantify GC content in 10 kb windows.
```
python slidingwindow_gc.py ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta 10000 10000 CroVir_genome.GC.10kb.txt
```

Run `slidingwindow_cpg.py` to quantify CpG content in 10 kb windows.
```
python slidingwindow_cpg.py ../genome_crotalus/CroVir_genome_L77pg_16Aug2017.fasta 10000 10000 CroVir_genome.CpG.10kb.txt
```

Compare distributions of GC and CpG content on autosomes and the sex chromosomes in `W_chromosome_composition.R`.

### 2. GC and CpG content in other amniotes

Here, compare autosomal and sex-linked GC content between prairie rattlesnake, birds, and mammals.

#### Set up environment
```
cd ./
mkdir genome_taenopygia
mkdir genome_homo
mkdir genome_mus
```

#### Retrieve genome data
```
cd ./genome_taeniopygia
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna.gz
gunzip GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna.gz 

cd ../genome_homo
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz 

cd ../genome_mus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
```

The chicken genome is already in the `./genome_gallus/` directory.

#### Calculate GC content for each genome in 10 kb windows
```
python slidingwindow_gc_content.py ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.fna 10000 10000 ./gc/genome_gallus.GC.10kb.txt
python slidingwindow_gc_content.py ./genome_taeniopygia/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna 10000 10000 ./gc/genome_taeniopygia.GC.10kb.txt
python slidingwindow_gc_content.py ./genome_homo/GCF_000001405.39_GRCh38.p13_genomic.fna 10000 10000 ./gc/genome_homo.GC.10kb.txt
python slidingwindow_gc_content.py ./genome_mus/GCF_000001635.27_GRCm39_genomic.fna 10000 10000 ./gc/genome_mus.GC.10kb.txt
```

#### Calculate CpG content for each genome in 10 kb windows
```
python slidingwindow_cpg.py ./genome_gallus/GCF_000002315.6_GRCg6a_genomic.fna 10000 10000 ./gc/genome_gallus.CpG.10kb.txt
python slidingwindow_cpg.py ./genome_taeniopygia/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna 10000 10000 ./gc/genome_taeniopygia.CpG.10kb.txt
python slidingwindow_cpg.py ./genome_homo/GCF_000001405.39_GRCh38.p13_genomic.fna 10000 10000 ./gc/genome_homo.CpG.10kb.txt
python slidingwindow_cpg.py ./genome_mus/GCF_000001635.27_GRCm39_genomic.fna 10000 10000 ./gc/genome_mus.CpG.10kb.txt
```

#### Analysis of GC/CpG content across amniotes

Run `GC_CpG_content_autosomes_sex-chromosomes.R` to compare distributions of GC and CpG content on autosomes and sex chromosomes across the sampled species.

## Repeat content

From the annotation of repeats it is clear that mdg4 elements contribute a substantial proportion of the total repetitive sequence on the W chromosome, which is over 80% repetitive. With this provisional result, we will examine variation in overall repeat content across the genome, along with a focus on mdg4 element abundance.

### Set up environment
```
mkdir repeats
```

### Quantify overall repeat and mdg4 content on W scaffolds

First, run `fasta_seq_length.py` on the W chromosome assembly to output a scaffold list with lengths.
```
python fasta_seq_length.py ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta > ./resources/Cviridis_CV0650_candidate_W.rescaffold.rename.lengths.txt
```

Then, run `scaffold_mdg4_repeat_content.py` to quantify repeat content per W scaffold.
```
python scaffold_mdg4_repeat_content.py ./resources/Cviridis_CV0650_candidate_W.rescaffold.rename.lengths.txt ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.gff3 ./repeats/Cviridis_CV0650_candidate_W.repeat_content_mdg4.txt
```

### Quantify overall mdg4 content on autosomes and the Z chromosome

As with GC content, these analyses will be performed in 10 kb sliding windows. First, extract mdg4 annotations from male reference annotation (these were annotated as their unfortunate synonym 'Gypsy' elements based on RepeatMasker libraries).
```
touch ./repeats/CroVir_genome_L77pg_16Aug2017.final.reformat.mdg4.sort.gff; grep "Gypsy" ./genome_crotalus/CroVir_genome_L77pg_16Aug2017.final.reformat.repeat.masked.sort.gff >> ./repeats/CroVir_genome_L77pg_16Aug2017.final.reformat.mdg4.sort.gff
```

Then run `window_mdg4_repeat_content.py` to calculate mdg4 content in 10 kb windows, using BED file for the male reference in `./resources/`.
```
python window_mdg4_repeat_content.py ./resources/CroVir_Dovetail_10Kb_window.bed ./repeats/CroVir_genome_L77pg_16Aug2017.final.reformat.mdg4.sort.gff ./repeats/CroVir_genome.repeat_content_mdg4.10kb.txt
```

## mdg4 retroelement GC-richness

The large abundance of mdg4 elements on the W chromosome may influence its overall nucleotide composition (i.e., GC-richness) if these elements also tend to be GC-rich.

### 1. GC content of repeat elements in snake library

First, examine the GC content of non-redundant elements in the snake transposable element library used to annotate the genome, which is in `./resources/Snakes_Known_TElib.fasta`.
```
python scaffold_gc.py ./resources/Snakes_Known_TElib.fasta ./repeats/Snakes_Known_TElib.GC.txt
```

Parse mdg4 results and remove redundant results.
```
touch ./repeats/Snakes_Known_TElib.mdg4.GC.txt; echo -e "chrom\tGC_content" >> ./repeats/Snakes_Known_TElib.mdg4.GC.txt; grep 'Gypsy' ./repeats/Snakes_Known_TElib.GC.txt | awk '!seen[$0] {print} {++seen[$0]}' >> Snakes_Known_TElib.mdg4.GC.txt
```

Replace '#' and '/' in output to make results readable in R.
```
sed -i 's/#/_/' ./repeats/Snakes_Known_TElib.GC.txt; sed -i.bak 's./._.g' ./repeats/Snakes_Known_TElib.GC.txt
sed -i 's/#/_/' ./repeats/Snakes_Known_TElib.mdg4.GC.txt; sed -i.bak 's./._.g' ./repeats/Snakes_Known_TElib.mdg4.GC.txt
```

Compare mdg4 element GC to all elements using `./R/mdg4_element_GC_content.R`.

### 2. GC-richness of annotated mdg4 elements on the W chromosome & correlations between GC content and mdg4 density

Convert repeat annotation GFF to BED format.
```
grep "^[^#]" ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.gff3 | awk 'BEGIN{OFS="\t"} {print $1,$4-1,$5,$9}' - > ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.bed
```

Use bedtools `getfasta` to extract sequence of each annotated repeat element. The -name flag specifies to name the output sequence based on the fourth input column.
```
bedtools getfasta -fi ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta -bed ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.bed -name -fo ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.fasta
```

Run `scaffold_gc.py` to calculate GC content per extracted sequence.
```
python scaffold_gc.py ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.fasta ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt
```

Extract entries for mdg4, L1, and L2 elements.
```
touch ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_mdg4.GC.txt; echo -e 'chrom\tGC_content' >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_mdg4.GC.txt; grep 'Gypsy' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_mdg4.GC.txt
touch ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L1-CIN4.GC.txt; echo -e 'chrom\tGC_content' >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L1-CIN4.GC.txt; grep 'L1' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L1-CIN4.GC.txt; grep 'CIN4' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L1-CIN4.GC.txt
touch ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L2-CR1-Rex.GC.txt; echo -e 'chrom\tGC_content' >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L2-CR1-Rex.GC.txt; grep 'L2' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L2-CR1-Rex.GC.txt; grep 'CR1' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L2-CR1-Rex.GC.txt; grep 'Rex' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt >> ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L2-CR1-Rex.GC.txt
```

Filter all repeats output to remove simple sequence repeats.
```
grep -v ')n' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt > ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.filter.GC.txt
```

Remove mdg4 elements from background for comparison with non-mdg4 elements.
```
grep -v 'Gypsy' ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt > ./repeats/Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_no-mdg4.GC.txt
```

Run `./R/mdg4_element_GC_content.R` to compare distributions of GC content and to calculate correlation coefficients between mdg4 density and GC content across W-linked scaffolds.

## Annotation of full-length LTR elements

### 1. Identification of full-length LTR retroelements on the W chromosome

Comparisons of full-length retroelements across the genome will be useful for understanding if the W chromosome acts as a 'refugium' for active or recently active elements capable of self-replication. Identification of full-length elements will make use of LTRharvest, which is part of the GenomeTools suite.

#### Set up environment
```
mkdir ./repeats/LTRharvest
mkdir ./repeats/LTRharvest/hmm
mkdir ./repeats/LTRharvest/fasta
cd ./repeats/LTRharvest
```

#### Retrieve and format search databases

The databases are available here:
* [Pfam](http://pfam.xfam.org/) - downloads are available under the `FTP` tab -> `current_release`.
* [GyDB](https://gydb.org/index.php/Collection_HMM)

Move hmm files from databases into `hmm` subdirectory.
```
wget https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
unzip GyDB_collection.zip
cd GyDB_collection/profiles
cp *.hmm ../../
cd ../../
rm GyDB_collection.zip

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
```

#### Prepare input data and run LTRharvest

First, generate an index for the input sequence data.
```
gt suffixerator -db ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta -indexname Cviridis_CV0650_candidate_W.rescaffold.rename.fasta -tis -suf -lcp -des -ssp -sds -dna
```

Run LTRharvest to prepare the genome data.
```
gt ltrharvest -index ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta -gff3 Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.gff -out Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.ltr.fa
gt ltrharvest -index ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta -seqids yes -tabout no > Cviridis_CV0650_candidate_W.rescaffold.rename.ltrharvest.out
```

Run LTRdigest to extract full-length LTR elements.
```
gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -checkids yes Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.gff > Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.sorted.gff
mv Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.sorted.gff Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.gff
gt ltrdigest -hmms ./hmm/*.hmm -aaout yes -outfileprefix Cviridis_CV0650_candidate_W.rescaffold.rename.fasta_ltrdigest Cviridis_CV0650_candidate_W.rescaffold.rename.fasta.gff ./resources/annotation/Cviridis_CV0650_candidate_W.rescaffold.rename.fasta > Cviridis_CV0650_candidate_W.rescaffold.rename.fasta_ltrdigest_output_gff
```

Move individual digest sequence files to the `fasta` subdirectory.
```
mv *.fas ./fasta
```

Extract GFF entries for full-length LTRs.
```
grep -v '#' Cviridis_CV0650_candidate_W.rescaffold.rename.fasta_ltrdigest_output_gff | awk '$3 == "LTR_retrotransposon" {print}' > fl-LTR.Cviridis_CV0650_candidate_W.rescaffold.rename.gff
```

### 2. Identification of full-length LTRs on autosomes and the Z chromosome

#### Prepare input data and run LTRharvest

Generate index for sequence file.

```
gt suffixerator -db CroVir_genome_L77pg_16Aug2017.final_rename.fasta -indexname CroVir_genome_L77pg_16Aug2017.final_rename.fasta -tis -suf -lcp -des -ssp -sds -dna
```

Run LTRharvest to prepare the sequence file.
```
gt ltrharvest -index CroVir_genome_L77pg_16Aug2017.final_rename.fasta -gff3 CroVir_genome_L77pg_16Aug2017.final_rename.fasta.gff -out CroVir_genome_L77pg_16Aug2017.final_rename.fasta.ltr.fa
gt ltrharvest -index CroVir_genome_L77pg_16Aug2017.final_rename.fasta -seqids yes -tabout no > CroVir_genome_L77pg_16Aug2017.final_rename.ltrharvest.out
```

Run LTRdigest to extract full-length LTRs.
```
gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -checkids yes CroVir_genome_L77pg_16Aug2017.final_rename.fasta.gff > CroVir_genome_L77pg_16Aug2017.final_rename.fasta.sorted.gff
mv CroVir_genome_L77pg_16Aug2017.final_rename.fasta.sorted.gff CroVir_genome_L77pg_16Aug2017.final_rename.fasta.gff
gt ltrdigest -hmms ./hmm/*.hmm -aaout yes -outfileprefix CroVir_genome_L77pg_16Aug2017.final_rename.fasta_ltrdigest CroVir_genome_L77pg_16Aug2017.final_rename.fasta.gff CroVir_genome_L77pg_16Aug2017.final_rename.fasta > CroVir_genome_L77pg_16Aug2017.final_rename.fasta_ltrdigest_output_gff
```

Move individual digest sequence files to the `fasta` subdirectory.
```
mv *.fas ./fasta
```

Extract GFF entries for full-length LTRs.
```
grep -v '#' CroVir_genome_L77pg_16Aug2017.final_rename.fasta_ltrdigest_output_gff > tmp.output.gff
for i in $(seq 0 17); do grep -w seq${i} tmp.output.gff >> fl-LTR.CroVir_genome_L77pg_16Aug2017.ChromAssigned.gff; done
rm tmp.output.gff
```

Note: comparisons of full-length LTRs including these outputs will be restricted to the chromosome-assigned scaffolds.

## Refugium & toxicity index analyses

Calculation of the refugium and toxicity indices for all repeat elements and specific classes is done based on the annotations described in the sections above. Quantitative details from these annotations are found in the [online supplement](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/gbe/14/9/10.1093_gbe_evac116/2/evac116_supplementary_data.zip?Expires=1667927995&Signature=C78XU3xHv7GsI2uW83wfbVxvRVOtPcfZ-xjCNYv5SatclEk3N3Pxh6cOCu5eeULt732NxEIVDSv8ElkiQ7vqtFOLb5toHRLE755974IUKHH9ZhKH5QmUNyW9mqRKcd2pX7aHPJwN8EGgAvAppGvgYZ0IGxAm-QpNNwNp94eP2KhJtAHfIHsDagLKLP3~oGcAiZI9eFoj5tZH8YVSK4pD3IIROWfcohbGRS08oxoNTYeN1kdYbTO2R7DyJcQ6aTbjXFj4o7fQjyGFAIJ241O2OKLkb2TlKe8C2GZsxlgXWJYOHsea7LaGGNRAl0r28iNhN0JP2coejJ1mDGcnsz~BnQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) for the paper.

To perform the statistical tests of the refugium hypothesis, run `./R/repeat_refugium_hypothesis_tests.R`.

## Gene expression analysis

For inquiries about the following gene expression analyses, please contact blair.perry[at]wsu.edu

#### 1. Map RNA-seq data to updated genome annotation
Quality trim raw RNA-seq reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). This step was also used to give each sample a simpler alias (i.e., Kid_24h_CV3_Cvv19_CAGGCG -> kidneyF1).
```bash
trimmomatic PE Cvv29_TACAGC_L005_R1_001.fastq.gz Cvv29_TACAGC_L005_R2_001.fastq.gz testes_forward_paired.fq.gz testes_forward_unpaired.fq.gz testes_reverse_paired.fq.gz testes_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Cvv35_AAACAC_L005_R1_001.fastq.gz Cvv35_AAACAC_L005_R2_001.fastq.gz ovaries_forward_paired.fq.gz ovaries_forward_unpaired.fq.gz ovaries_reverse_paired.fq.gz ovaries_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Kid_24h_CV3_Cvv19_CAGGCG_L005_R1_001.fastq.gz Kid_24h_CV3_Cvv19_CAGGCG_L005_R2_001.fastq.gz kidneyF1_forward_paired.fq.gz kidneyF1_forward_unpaired.fq.gz kidneyF1_reverse_paired.fq.gz kidneyF1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Kid_24h_CV5_Cvv18_CACTCA_L005_R1_001.fastq.gz Kid_24h_CV5_Cvv18_CACTCA_L005_R2_001.fastq.gz kidneyM1_forward_paired.fq.gz kidneyM1_forward_unpaired.fq.gz kidneyM1_reverse_paired.fq.gz kidneyM1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Kid_24h_CV6_Cvv20_CATGGC_L005_R1_001.fastq.gz Kid_24h_CV6_Cvv20_CATGGC_L005_R2_001.fastq.gz kidneyM2_forward_paired.fq.gz kidneyM2_forward_unpaired.fq.gz kidneyM2_reverse_paired.fq.gz kidneyM2_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Kid_24h_CV8_Cvv17_CACGAT_L005_R1_001.fastq.gz Kid_24h_CV8_Cvv17_CACGAT_L005_R2_001.fastq.gz kidneyF2_forward_paired.fq.gz kidneyF2_forward_unpaired.fq.gz kidneyF2_reverse_paired.fq.gz kidneyF2_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE liver_24h_CV3_Cvv08_CGTACG_L005_R1_001.fastq.gz liver_24h_CV3_Cvv08_CGTACG_L005_R2_001.fastq.gz liverF1_forward_paired.fq.gz liverF1_forward_unpaired.fq.gz liverF1_reverse_paired.fq.gz liverF1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE liver_24h_CV5_Cvv07_GTTTCG_L005_R1_001.fastq.gz liver_24h_CV5_Cvv07_GTTTCG_L005_R2_001.fastq.gz liverM1_forward_paired.fq.gz liverM1_forward_unpaired.fq.gz liverM1_reverse_paired.fq.gz liverM1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Liver_24h_CV6_Cvv06_GTGGCC_L005_R1_001.fastq.gz liver_24h_CV6_Cvv06_GTGGCC_L005_R2_001.fastq.gz liverM2_forward_paired.fq.gz liverM2_forward_unpaired.fq.gz liverM2_reverse_paired.fq.gz liverM2_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
trimmomatic PE Liver_24h_CV8_Cvv05_GTGAAA_L005_R1_001.fastq.gz liver_24h_CV8_Cvv05_GTGAAA_L005_R2_001.fastq.gz liverF2_forward_paired.fq.gz liverF2_forward_unpaired.fq.gz liverF2_reverse_paired.fq.gz liverF2_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
```

Convert updated GFF file (including original + new W chromosome annotation) to GTF format with [gffread](https://github.com/gpertea/gffread).
```bash
gffread CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.withWannot.sort.gff -T -F -o crovir_wWannot.gtf
```

Generate genome index files needed for mapping with [STAR](https://github.com/alexdobin/STAR).
```bash
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir ./reference/index --genomeFastaFiles ./reference/CroVir_genome_L77pg_16Aug2017.final_rename.withWscaffs.fasta --sjdbGTFfile ./reference/crovir_wWannot.gtf 
```

Map trimmed RNA-seq reads to the genome using STAR. Note, an example command is shown below for a female kidney sample; this command was run separately for each of the 10 RNA-seq samples used in these analyses. 
```bash
# Mapping with STAR (Example command for Kidney F1 sample)
STAR --genomeDir ./reference/index --runThreadN 12 --outFilterMultimapNmax 1 --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.1 --twopassMode Basic --sjdbGTFfile ./reference/ref1_origAndW/crovir_wWannot.gtf --quantMode GeneCounts --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR_mapped/kidneyF1 --readFilesIn ./trimmed/paired/kidneyF1_forward_paired.fq.gz ./trimmed/paired/kidneyF1_reverse_paired.fq.gz
```

Quantify reads using [featureCounts](http://subread.sourceforge.net/). Note that the GTF file generated by gffread above has gene identifiers in the 'transcript_id' field, and therefore the flag `-g transcript_id` in the following command results in gene-level counts (not transcript-level). 
```bash
featureCounts -p -t exon -g transcript_id -a ./reference/crovir_wWannot.gtf -o .raw_counts/cvv_ref1_rawCounts_12.18.21.txt STAR_mapped/*Aligned.sortedByCoord.out.bam 
```

#### 2. Assess "detectable" expression in male and female samples.
Assess whether genes are detected (raw count > 0) in male and female samples.
Note: The following code is all run within R. 

Load packages.
```R
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(IHW)
library(scico)
library(WebGestaltR)
```

Read in table with ZW gametolog information.
```R
zw.geneInfo <- readxl::read_xlsx('CVV_Wchrom_FunctionalAnalyses_March2022/W_genes_Zpositions.xlsx') %>% janitor::clean_names()
```

Read in raw counts from featureCounts and reformat into simple dataframe.
```R
autoZW.unique.raw_counts <- read_tsv('raw_counts/cvv_ref1_rawCounts_12.18.21.txt',comment = '#') %>% 
  janitor::clean_names() %>% # clean up header names
  select(-2,-3,-4,-5) # remove unnecessary columns

names(autoZW.unique.raw_counts) <- str_remove_all(names(autoZW.unique.raw_counts),'ref1_unique_star_mapped_|aligned_sorted_by_coord_out_bam|_aligned_sorted_by_coord_out_bam') # further clean up headers

autoZW.unique.raw_counts.simple <- as.data.frame(autoZW.unique.raw_counts[,-c(1,2)])
row.names(autoZW.unique.raw_counts.simple) <- autoZW.unique.raw_counts$geneid
```

Filter table to only genes on W chromosome, and tally the number of genes with "detected" expression (raw count > 0) only in females. 
```R
autoZW.unique.raw_counts.simple.Wonly <- as_tibble(autoZW.unique.raw_counts.simple,rownames = 'id') %>% 
  filter(str_detect(id,'scaffold-W')) %>% 
  filter(str_detect(id,'scaffold-W393371|trna',negate=T)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  mutate(total_female = rowSums(.[,c(2,3,6,7,10)])) %>% 
  mutate(total_male = rowSums(.[,c(4,5,8,9,11)]))

# Tally and plot number of W genes with detected expression in females
autoZW.unique.raw_counts.simple.Wonly %>% 
  group_by(total_female > 0) %>% 
  tally() %>% 
  ggplot(aes(x=`total_female > 0`,y=n)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=n),nudge_y = -4,color='white') +
  theme_linedraw()

# 137 expressed, 82 not
# 82 genes on the W with NO female sample expression counts 

# Tally number of W genes with detected expression ONLY in females
autoZW.unique.raw_counts.simple.Wonly %>% 
  filter(total_male == 0 & total_female > 0) %>% 
  group_by(total_female > 0) %>% 
  tally()

# 103 genes detected ONLY in female samples

# Tally number of W chrom genes with detected expression on each stratum
autoZW.unique.raw_counts.simple.Wonly %>% 
  mutate(id = str_remove_all(id,'-mRNA-1')) %>% 
  left_join(zw.geneInfo,by=c('id'='feature_abbrev')) %>% 
  filter(!is.na(stratum)) %>% 
  group_by(total_female > 0,stratum) %>% 
  tally()

# Older: 44 expressed, 19 not
# Recent: 45 expressed, 18 not
# Unknown: 48 expressed, 45 not
```

#### 3. Perform pairwise comparisons between male and female samples
Use [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to perform pariwise comparisons between male and female samples. 
```R
# Filter out genes with no expression across all samples
autoZW.unique.raw_counts.simple.filtered <- autoZW.unique.raw_counts.simple[rowSums( autoZW.unique.raw_counts.simple != 0 ) > 0,]

# Define male/female groups for pairwise comparisons
autoZW.unique.colData <- (DataFrame(condition=group <- factor(c('F','F','M','M','F','F','M','M','F','M'))))

# Set up DESeq2 dataset object
autoZW.unique.dds <- DESeqDataSetFromMatrix(autoZW.unique.raw_counts.simple.filtered,autoZW.unique.colData,formula(~condition))

# Perform DESeq2 analysis
autoZW.unique.dds <- DESeq(autoZW.unique.dds)
```

Extract pariwise comparisons results and perform [IHW](https://bioconductor.org/packages/release/bioc/html/IHW.html) p-value correction.
```R
# Generate data frame with pairwise comparison results
autoZW.unique.deResFemVsMale <- as.data.frame(results(autoZW.unique.dds, contrast=c('condition','F','M')))

# Calculate IHW p-value using baseMean as the covariate
autoZW.unique.ihwRes <- ihw(pvalue ~ baseMean,  data = autoZW.unique.deResFemVsMale, alpha = 0.05)

# Add IHW p-value to result table
autoZW.unique.deResFemVsMale$IHW_pvalue <- autoZW.unique.ihwRes@df$adj_pvalue
autoZW.unique.deResFemVsMale <- autoZW.unique.deResFemVsMale[order(autoZW.unique.deResFemVsMale$IHW_pvalue),]

# Make table of genes with female biased expression (significantly diff expressed with IHW p-value < 0.05 and log2FC > 0)
autoZW.unique.femUpreg <- as_tibble(autoZW.unique.deResFemVsMale,rownames = 'gene') %>% 
  filter(IHW_pvalue < 0.05 & log2FoldChange > 0)

nrow(autoZW.unique.femUpreg) # 33 significantly upregulated genes in female

# Tally number of female-biased genes that are located on the W chrom
autoZW.unique.femUpreg %>% 
  mutate(onW = ifelse(str_detect(gene,'scaffold-W'),T,F)) %>% 
  group_by(onW) %>% 
  tally()

# 31 out of 33 are on the W
```

Generate TPM normalized counts for plotting. 
```R
# Read in and clean raw count table
raw.counts <- read_tsv('raw_counts/raw_counts_ref1/cvv_ref1_rawCounts_12.18.21.txt',comment = '#') %>% 
  janitor::clean_names() %>% 
  select(-2,-3,-4,-5)

# Clean headers
names(raw.counts) <- str_remove_all(names(raw.counts),'ref1_unique_star_mapped_|aligned_sorted_by_coord_out_bam|_aligned_sorted_by_coord_out_bam')

# Get gene lengths from table
gene.lengths <- raw.counts$length

# Make matrix of raw counts
raw.counts.matrix <- as.matrix(raw.counts[,-c(1,2)])
row.names(raw.counts.matrix) <- raw.counts$geneid

# Calculate Transcripts per million (TPM) counts
count.over.length  <-  raw.counts.matrix / gene.lengths

tpm.counts <- as.data.frame(t( t(count.over.length) * 1e6 / colSums(count.over.length) )) %>% 
  rownames_to_column('id')  

# Format a count table that includes male and female mean, whether a gene is on the W chromosome, and if so, which stratum its located within.
tpm.counts.final <- as_tibble(tpm.counts) %>% 
  mutate(femaleMean = rowMeans(.[,c(2,3,6,7,10)]),
         maleMean = rowMeans(.[,c(4,5,8,9,11)])) %>% 
  mutate(id = str_remove_all(id,'-mRNA-1')) %>% 
  left_join(zw.geneInfo,by=c('id'='feature_abbrev')) %>% 
  arrange(stratum,z_chromosome_position) %>% 
  # select(-c(15,16)) %>% 
  select(id,gene, everything())

# Write TPM counts to file
write_tsv(tpm.counts.final,'_autoZW.allGene.TPMCounts_03.02.22.tsv')
```

#### 4. Plot heatmap of ZW gametolog expression
Read in TPM counts and filter to W chromosome genes only, add column indicating whether a gene is upregulated/biased in females.
```R
tpm.counts <- read_tsv('_autoZW.allGene.TPMCounts_03.02.22.tsv') %>% 
  filter(str_detect(id,'scaffold-W')) %>% 
  filter(!is.na(stratum)) %>% 
  mutate(stratum = factor(stratum,levels=c('recent','older','unknown'))) %>% 
  arrange(stratum,z_chromosome_position) %>% 
  mutate(female_upreg = ifelse(id %in% female.upreg$gene,1,0))
```

Make a dataframe of TPM counts and an annotation table indicating which stratum each gene is located within for use in heatmap plotting. 
```R
exp.zw.maleFemaleHeat <- tpm.counts %>% select(id,ovaries,contains('kidney_f'),contains('liver_f'),femaleMean,maleMean,testes,contains('kidney_m'),contains('liver_m')) %>% 
  as.data.frame() %>% 
  column_to_rownames('id') 
  
exp.zw.femaleHeat.annot <- tpm.counts %>% select(id,stratum,female_upreg) %>%  as.data.frame() %>% column_to_rownames('id')

# Filter out genes where the stratum is not known
exp.zw.femaleHeat.annot.noUnknown <- exp.zw.femaleHeat.annot %>% 
  filter(stratum != 'unknown')

exp.zw.maleFemaleHeat.noUnknown <- exp.zw.maleFemaleHeat %>% 
  filter(row.names(exp.zw.maleFemaleHeat) %in% 
```

Plot heatmap with [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html).
```R
pheatmap(log10(exp.zw.maleFemaleHeat.noUnknown+1),
         border_color = NA,
         scale='none',
         cellheight = 4,cellwidth = 10,
         cluster_rows = F,cluster_cols = F,
         color = scico(50,palette = 'lajolla',direction = -1),
         show_rownames = F,
         gaps_col = c(5,7),
         annotation_row = exp.zw.femaleHeat.annot.noUnknown)
```

#### 5. Characterize genes with female-biased expression using GO term analyses.

Load in foreground (genes on W chromsome) and background (genes not on W chromosome) gene lists.
```R
w.genes <- read_tsv('~/Downloads/Wchrom_GOtermAnalysis/WChromGenes.txt',col_names = 'symbol')
bg.genes <- read_tsv('~/Downloads/Wchrom_GOtermAnalysis/bg_forWChromGO.txt',col_names = 'symbol')
```

Perform GO and pathway overrepresentation analysis of W genes versus non-W genes using [Webgestalt](https://cran.r-project.org/web/packages/WebGestaltR/index.html).
```R

go.res <- WebGestaltR(enrichMethod = 'ORA',
            enrichDatabase = c('geneontology_Biological_Process_noRedundant','geneontology_Molecular_Function_noRedundant','geneontology_Cellular_Component_noRedundant'),
            interestGene = w.genes$symbol,
            interestGeneType = 'genesymbol',
            referenceGene = bg.genes$symbol,
            referenceGeneType = 'genesymbol',
            sigMethod = 'top',
            topThr = 20,
            isOutput = F)

go.res$analysis <- 'Gene Ontology'

path.res <- WebGestaltR(enrichMethod = 'ORA',
             enrichDatabase = c('pathway_KEGG','pathway_Panther','pathway_Reactome','pathway_Wikipathway'),
             interestGene = w.genes$symbol,
             interestGeneType = 'genesymbol',
             referenceGene = bg.genes$symbol,
             referenceGeneType = 'genesymbol',
             sigMethod = 'top',
             topThr = 20,
             isOutput = F)

path.res$analysis <- 'Pathway'

```

Combine and plot results.
```R
all.res <- go.res %>% 
  bind_rows(path.res) %>% 
  arrange(analysis,overlap) %>% 
  mutate(description = factor(description,levels=unique(.$description))) %>% 
  mutate(signif = ifelse(FDR < 0.05,T,F)) %>% 
  mutate(database = str_split_fixed(database,'[_]',2)[,2] %>% str_replace_all('_',' ') %>% str_to_title())

ggplot(all.res,aes(x=overlap,y=description,alpha=signif,fill=database)) +
  geom_bar(stat='identity',orientation = 'y') + geom_point(data=subset(all.res,signif==TRUE),aes(x=overlap*1.1,y=description),pch='*',size=8,show.legend = F) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(all.res$overlap)+2)) +
  scale_alpha_manual(values = c('TRUE'=1,'FALSE'=0.35),guide='none')+
  facet_wrap(~analysis,ncol = 1,scales='free_y') +
  xlab('Number of Genes') +
  ylab('Term') +
  labs(alpha='FDR < 0.05',fill='Database') +
  theme_classic() + theme(strip.background = element_rect(fill = 'grey30',colour = 'NA'),strip.text = element_text(color='white'))
```

#### 6. GO term overrepresentation analysis of duplicated and translocated genes.
The following analyses are performed in R. 

Set up environment and read in data.
```R
library(tidyverse)
library(WebGestaltR)

w.dups.trans <- readxl::read_xlsx('W_duplication_translocation_genes.xlsx') %>% 
  janitor::clean_names()

w.dupsList <- w.dups.trans %>% 
  select(w_specific_duplicated_genes) %>% 
  filter(!is.na(w_specific_duplicated_genes))

w.transList <- w.dups.trans %>% 
  select(translocated_genes) %>% 
  filter(!is.na(translocated_genes))

bg.genes <- read_tsv('~/Downloads/Wchrom_GOtermAnalysis/bg_forWChromGO.txt',col_names = 'symbol')
```

Perform GO analysis for duplicated genes using [WebgestaltR](https://cran.r-project.org/web/packages/WebGestaltR/index.html).
```R
dups.go.res <- WebGestaltR(enrichMethod = 'ORA',
                           enrichDatabase = c('geneontology_Biological_Process_noRedundant','geneontology_Molecular_Function_noRedundant','geneontology_Cellular_Component_noRedundant'),
                           interestGene = w.dupsList$w_specific_duplicated_genes,
                           interestGeneType = 'genesymbol',
                           referenceGene = bg.genes$symbol,
                           referenceGeneType = 'genesymbol',
                           sigMethod = 'top',
                           topThr = 20,
                           isOutput = F)

dups.go.res$analysis <- 'Gene Ontology'

dups.path.res <- WebGestaltR(enrichMethod = 'ORA',
                        enrichDatabase = c('pathway_KEGG','pathway_Panther','pathway_Reactome','pathway_Wikipathway'),
                        interestGene = w.dupsList$w_specific_duplicated_genes,
                        interestGeneType = 'genesymbol',
                        referenceGene = bg.genes$symbol,
                        referenceGeneType = 'genesymbol',
                        sigMethod = 'top',
                        topThr = 20,
                        isOutput = F)

dups.path.res$analysis <- 'Pathway'

dups.all.res <- dups.go.res %>% 
  bind_rows(dups.path.res) %>% 
  arrange(overlap) %>% 
  mutate(description = factor(description,levels=unique(.$description))) %>% 
  mutate(signif = ifelse(FDR < 0.05,T,F)) %>% 
  mutate(database = str_split_fixed(database,'[_]',2)[,2] %>% str_replace_all('_',' ') %>% str_to_title())

ggplot(dups.all.res,aes(x=overlap,y=description,alpha=signif,fill=database)) +
  geom_bar(stat='identity',orientation = 'y') +
  geom_point(data=subset(dups.all.res,signif==TRUE),aes(x=overlap*1.1,y=description),pch='*',size=8,show.legend = F) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(dups.all.res$overlap)+2)) +
  scale_alpha_manual(values = c('TRUE'=1,'FALSE'=0.5),guide='none')+
  facet_wrap(~analysis,ncol = 1,scales='free_y') +
  xlab('Number of Genes') +
  ylab('Term') +
  labs(alpha='FDR < 0.05',fill='Database') +
  theme_classic() + theme(strip.background = element_rect(fill = 'grey30',colour = 'NA'),strip.text = element_text(color='white'))
```

Perform GO analysis for translocated genes using [WebgestaltR](https://cran.r-project.org/web/packages/WebGestaltR/index.html).
```R
trans.go.res <- WebGestaltR(enrichMethod = 'ORA',
                           enrichDatabase = c('geneontology_Biological_Process_noRedundant','geneontology_Molecular_Function_noRedundant','geneontology_Cellular_Component_noRedundant'),
                           interestGene = w.transList$translocated_genes,
                           interestGeneType = 'genesymbol',
                           referenceGene = bg.genes$symbol,
                           referenceGeneType = 'genesymbol',
                           sigMethod = 'top',
                           topThr = 20,
                           isOutput = F)

trans.go.res$analysis <- 'Gene Ontology'

trans.path.res <- WebGestaltR(enrichMethod = 'ORA',
                             enrichDatabase = c('pathway_KEGG','pathway_Panther','pathway_Reactome','pathway_Wikipathway'),
                             interestGene = w.transList$translocated_genes,
                             interestGeneType = 'genesymbol',
                             referenceGene = bg.genes$symbol,
                             referenceGeneType = 'genesymbol',
                             sigMethod = 'top',
                             topThr = 20,
                             isOutput = F)

trans.path.res$analysis <- 'Pathway'


trans.all.res <- trans.go.res %>% 
  bind_rows(trans.path.res) %>% 
  arrange(overlap) %>% 
  mutate(description = factor(description,levels=unique(.$description))) %>% 
  mutate(signif = ifelse(FDR < 0.05,T,F)) %>% 
  mutate(database = str_split_fixed(database,'[_]',2)[,2] %>% str_replace_all('_',' ') %>% str_to_title())

ggplot(trans.all.res,aes(x=overlap,y=description,alpha=signif,fill=database)) +
  geom_bar(stat='identity',orientation = 'y') +
  geom_point(data=subset(trans.all.res,signif==TRUE),aes(x=overlap*1.1,y=description),pch='*',size=8,show.legend = F) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(trans.all.res$overlap)+2)) +
  scale_alpha_manual(values = c('TRUE'=1,'FALSE'=0.5),guide='none')+
  facet_wrap(~analysis,ncol = 1,scales='free_y') +
  xlab('Number of Genes') +
  ylab('Term') +
  labs(alpha='FDR < 0.05',fill='Database') +
  theme_classic() + theme(strip.background = element_rect(fill = 'grey30',colour = 'NA'),strip.text = element_text(color='white'))
```

## W-specific gene duplications

UNDER CONSTRUCTION

## Sex-linked divergence between pitvipers

UNDER CONSTRUCTION

## ZW gametolog GC3 analysis

UNDER CONSTRUCTION

## Appendix: Indian cobra analysis

UNDER CONSTRUCTION










































