---
title: Read alignment
type: docs
prev: docs/folder/2-prevariant-filtering
next: docs/folder/4-variant-calling
---

## Aligning reads to a genome
Once the reads look good, then we need to align them to a genome so that we can make sense of them. For this, we need (1) a program to use, (2) the reads we want to align, and (3) the genome that we want to align the reads to. The common alignment programs that are used for Illumina short reads include [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) and [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). These are both available as modules: `module load bwa` // `module load bowtie2`.

With your own data in the future, you may find that other aligners may be better for your purposes, depending on what kind of data you’re working with. There also are “workflow” programs that can be used to do many of these steps together: [STACKS](https://catchenlab.life.illinois.edu/stacks/manual/), [ipyrad](https://ipyrad.readthedocs.io/en/master/index.html), etc. (but these usually use bwa or bowtie2 under the hood). These can be nice because they don't require as much hands-on work, but they can also sometimes become "black boxes" that can get you into trouble if you're not thinking explicitly about the settings that you're using.

### Preparing the reference genome
We also need a reference genome to use! Look on NCBI genome database if you don’t have one, you can download the genome file (for our in-class work, we’ll be using the Arctic charr reference genome: [https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002910315.2/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002910315.2/), which I’ve already downloaded to our shared directory). I’ve already indexed the genome, but to do so with your own genome, you’ll need to use `bwa index genome.fa.gz` and `bowtie2-build genome.fa.gz genome_prefix` (these can take a while, depending on the size of your genome). Indices are necessary to allow the programs to more quickly access the genome once they start trying to align reads.

### Aligning single-end reads
For our work in class, we'll be working with the reads that have been assigned to a single individual; in practice, this process will have to be repeated separately with all of the individuals that you are working with. Let's say that we have our trimmed and filtered reads for an individual in the file `ind1_reads.fastq` and the reference genome that we're using is in `/xdisk/jrick/consbioinf/shared_data/Salvelinus_spp_genome.fasta`. We can then align these reads to the reference using:

```sh
bwa mem /xdisk/jrick/consbioinf/shared_data/Salvelinus_spp_genome.fasta ind1_reads.fastq > aln-ind1-bwa.sam
```

Or, if we'd prefer using bowtie2, then we can do the alignment using:

```sh
bowtie2 -x /xdisk/jrick/consbioinf/shared_data/Salvelinus_spp_genome -U ind1_reads.fastq -S aln-ind1-bowtie.sam
```
These processes both produce `.sam` files, which stands for **S**equence **A**lignment **M**ap file. These are human-readable files that contain information about where and how well the reads align to the genome. These are also generally *very large* files, and so we usually want to then convert them to `.bam` files, or **B**inary **A**lignment **M**ap files, which is a format in which the information is compressed. To do so, we'll use the program [SAMtools](https://www.htslib.org/doc/samtools.html), which is also installed as a module:

```sh
module load samtools

samtools view --bam --with-header --output aln-ind1-bwa.bam aln-ind1-bwa.sam
```
From here, many downstream applications will require that your `.bam` file is sorted, such that the reads are in the order that they occur in the genome. This can be done using the following code:

```sh
samtools sort --output aln-ind1-bwa.sorted.bam aln-ind1-bwa.bam
```
Once we have our sorted bamfile, then we are ready for the next step of variant calling!

### Aligning paired-end reads
A "paired-end" or "mate-pair" read consists of pair of mates, called mate 1 and mate 2. Pairs come with a prior expectation about (a) the relative orientation of the mates, and (b) the distance separating them on the original DNA molecule. Exactly what expectations hold for a given dataset depends on the lab procedures used to generate the data. For example, a common lab procedure for producing pairs is Illumina's Paired-end Sequencing Assay, which yields pairs with a relative orientation of FR ("forward, reverse"). This protocol usually yields pairs where the expected genomic distance from end to end is about 200-500 base pairs, and in some cases the reads will overlap in the middle.

When working with paired-end reads, the program calls are similar, but you will need to provide one "forward read" fastq file and one "reverse read" fastq file. You will also need to make a decision about whether you want all reads kept, even if they do not align to the same location as their partner, or if you want to discard any reads that align to a separate place in the genome.

```sh
bwa mem /xdisk/jrick/consbioinf/shared_data/Salvelinus_spp_genome.fasta ind1_reads_fwd.fastq ind1_reads_rev.fastq > aln-ind1-paired-bwa.sam
```
