# Rattlesnake W chromosome analysis

This repository contains details on the data processing and analysis steps used to assemble, annotate, and characterize the female-specific W chromosome, including comparative genomic analyses with other caenophidian snake and amniote species. This workflow is a companion to the methods described in Schield et al. (in review).

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
* [MultiQC](xxxx)
* [10x Genomics Supernova](xxxx) (v2.1.1)
* [NCBI BLAST](xxxx)
* [MashMap](xxxx)
* [Trinity](xxxx)
* [Maker](xxxx)
* [Agouti](xxxx)
* [Repeatmasker](xxxx)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bwa](http://bio-bwa.sourceforge.net/)
* [htslib](http://www.htslib.org/)
* [samtools](http://www.htslib.org/)
* [bgzip](http://www.htslib.org/)
* [tabix](http://www.htslib.org/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [mosdepth](xxxx)
* [STAR](xxxx)
* [R](https://cran.r-project.org/)

Note, I installed a number of these programs to my [conda](https://docs.conda.io/en/latest/) environment.

## Female genome assembly

### 10x Genomics Chromium linked-read sequencing

We generated a 10x Genomics Chromium library for a female prairie rattlesnake from snap-frozen liver tissue. The library was sequenced on an Illumina NovaSeq 6000 using 150 bp paired-end reads.

### FastQC/MultiQC analysis

Assess read quality using FastQC.

#### Set up environment







### Assembly

## Identification of W chromosome scaffolds

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

