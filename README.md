# Rattlesnake W chromosome analysis

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

1,500,000,000 x 56 / 150 = 560,000,000 input reads.

```
supernova run --maxreads 560000000 --id=CV0650_female_Cviridis --fastqs=./fastq/
```

Generate pseudohaplotype fasta assembly with a minimum scaffold length of 10 kb.

```
supernova mkoutput --style=pseudohap2 --asmdir=./CV0650_female_Cviridis/outs/assembly/ --outprefix=CV0650_female_Cviridis --minsize=10000 --headers=short 
```

## Identification of W chromosome scaffolds

Use homology with the male prairie rattlesnake reference genome and relative female/male read depths on female scaffolds to identify candidate W chromosome scaffolds.

### Set up environment

```
mkdir genome_crotalus
mkdir genome_crotalus_female
```


You can find the male 



## W chromosome annotation

## ZW gametolog divergence

### Identification of 1:1 ZW gametologs

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

