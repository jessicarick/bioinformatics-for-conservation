---
title: 14 - Molecular Forensics
type: docs
prev: docs/folder/13-hybridization
next: docs/folder/14-molecular-forensics
---

## Molecular Forensics
Molecular forensics is a generic term that is used to describe analyses of samples that have an unknown identity, where the goal is to determine the identity of the sample(s). 

*note that this page is a work in progress and currently incomplete*

There are many programs and pipelines out there for analyzing metabarcoding and metagenomic data. Some programs or pipelines that are commonly used with these types of data:
* [QIIME 2](https://qiime2.org/)
* [DADA2 R package](https://benjjneb.github.io/dada2/)
* [OBItools](https://git.metabarcoding.org/obitools/obitools/-/wikis/home)
* [metabaR R package](https://metabarfactory.github.io/metabaR/index.html)
* [USEARCH](https://www.drive5.com/usearch/)
* [VSEARCH](https://github.com/torognes/vsearch)

## eDNA Metabarcoding Example
For our in-class exercise, we are going to be analyzing an eDNA sample that was taken from an agave flower with the purpose of detecing *Leptonycteris* bat pollinators as a part of the study published in [Walker et al. 2022](https://doi.org/10.3390/ani12223075). Instead of specifically looking for bat DNA in the sample, we are going to analyze the data to see what all taxa we can find. For today, we'll be working completely on the UA HPC.

```sh
## usearch isn't a module, so we'll store the path to where it is installed
alias usearch='/xdisk/jrick/programs/usearch11.0.667_i86linux32'

## first, join paired reads
## note that in class, we had poor luck with the paired reads,
## so we decided to just use the R1 reads instead
usearch -fastq_join /xdisk/jrick/consbioinf/shared_data/week13_data/SRR21894625_R1.fastq \
	-reverse /xdisk/jrick/consbioinf/shared_data/week13_data/SRR21894625_R2.fastq \
	-fastqout SRR21894625_join.fastq

## filter specifically for our analyses
## to remove low quality sequences
usearch -fastq_filter /xdisk/jrick/consbioinf/shared_data/week13_data/SRR21894625_R1.fastq \
	-fastq_maxee 1.0 \ # can change to make more/less stringent
	-relabel Filt \
	-fastaout SRR21894625_R1_filtered.fa

## Find unique read sequences and their abundances
usearch -fastx_uniques SRR21894625_R1_filtered.fa \
	-relabel Uniq \
	-sizeout \
	-fastaout SRR21894625_R1_filtered_uniq.fa

## Make 97% similarity OTUs and filter chimeras
usearch -cluster_otus SRR21894625_R1_filtered_uniq.fa \
	-otus SSRR21894625_R1_filtered_otus.fa \
	-relabel Otu

## Denoise
usearch -unoise3 SRR21894625_R1_filtered_uniq.fa \
	-zotus SRR21894625_R1_zotus.fa

## Search against BOLD COI database
usearch -usearch_global SRR21894625_R1_zotus.fa \
	-db /xdisk/jrick/consbioinf/shared_data/week13_data/bold_COI_ref.udb  \
	-id 0.6 \
	-alnout SRR21894625_R1_bold_COI_hits.aln \
	-strand plus \
	-maxaccepts 3

## inspect the taxa in our results (to exit the "less" command, type "q")
less SRR21894625_bold_COI_hits.aln

```
