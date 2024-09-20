---
title: 2 - Quality Assessment and Filtering
type: docs
prev: docs/folder/1-ua-hpc-intro
next: docs/folder/3-read-alignment
---

## Why do we need to filter our data?

Sequencing is not a perfect process, and errors/biases can be introduced at many steps. As we read about in [Hemstrom et al. (2024)](https://www.nature.com/articles/s41576-024-00738-6), there are errors that get introduced through various steps in the sequencing process, including biases introduced during sequencing library preparation, polymerase errors during replication, errors that occur during the sequencing process itself, and inaccurate base calling. In addition, errors can be introduced when sequences are aligned to a reference genome.

“Filtering is an issue of paramount importance because every genomic dataset must be filtered, often repeatedly, and the same dataset filtered in different ways can yield entirely different results” (Hemstrom et al. 2024)

Additional reading on the subject: Shafer et al. https://doi.org/10.1111/2041-210X.12700 and Nazareno & Knowles https://doi.org/10.3389/fpls.2021.677009 


## How do we assess data quality?
When we get our hands on sequencing data, we always want to start with quality assessment, followed by read trimming and filtering, depending on how those data look. The most common program used to do quality assessment of raw `fastq` files is [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). This program can be used to get an idea of how the reads look and what kind of filtering may be needed. On the UA HPC, Fastqc is installed as a module, so we can use it as follows:

```sh
module load fastqc

fastqc -o output_directory/ sequencing_data.fastq
```

This will take a while to run on a whole lane of data, so for in-class purposes we will be working with a small subset of a seuqencing lane, `small_AC_1.clean.fastq`, which can be found in our shared `/xdisk/jrick/consbioinf/shared_data` directory. To run Fastqc on these data, we will use:

```sh
mkdir fastqc_output/
fastqc -o fastqc_output/ small_AC_1.clean.fastq
```

Once it runs, we can download the html and view the results. The easiest way to do this is to navigate to [https://ood.hpc.arizona.edu/](https://ood.hpc.arizona.edu/), navigate to the File Viewer, and then navigate to the `/xdisk/jrick/shared_data/fastqc_output/` directory. There will be a `html` file with the `small_AC_1.clean` prefix, and you can download this and then open it on your computer. 

As we look through the results, the program will flag what it thinks may be concerning to you; however, what things are actually concerning to you will depend on knowing your expectations for your data---for example, I know that the first 8-10 bases in our reads are the barcodes, so it doesn’t worry me that the per-base sequence content is doing funky things around there. However, we do find that there is an overrepresentation of Illlumina adapters, which suggests that we should probably trim these from our reads (but see note below that you may not actually need to trim these).

## How do we trim and filter our reads?
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) is commonly used for filtering/trimming sequencing reads
- Load on the UA HPC using `module load cutadapt`
- We probably want to quality trim the end of the reads (`-q 10` trims anything with quality < 10) and maybe trim the adapter sequences (if needed, using `-a ADAPTER` or `-a file:adapters.fasta`). In class, we chose to trim any bases with a quality < 20. To do this, we used the following code:

```sh
cutadapt -q 20 -o AC1-AC16_0505_009-q20.fastq /xdisk/jrick/consbioinf/shared_data/char_fastq/AC1-AC16_0505_009.fastq
```
  
- Info from cutadapt documentation: Under some circumstances, you may want to consider not trimming adapters at all. For example, a good library prepared for exome, genome or transcriptome sequencing should contain very few reads with adapters anyway. Also, some read mapping programs including BWA-MEM and STAR will soft-clip bases at the 3’ ends of reads that do not match the reference, which will take care of adapters implicitly.
- Cutadapt can also be used to demultiplex when many individuals are included in the same file (some sequencing facilities will give you demultiplexed fastq files; others, you’ll have to demultiplex yourself) — `cutadapt -e 1 -g ^file:AC_1_barcode_key.fasta -o “AC1_demux/AC1-demux-{name}.fastq” AC_1.clean.fastq`

[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) is an alternative to cutadapt that can be used to trim reads (removing adapters, removing poor quality bases at the ends of reads, etc) and it has lots of options!
- Load using `module load trimmomatic`
- Example code: `trimmomatic SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 MINLEN:36` -- where `2:30:10` refers to the number of seed mismatches (the maximum mismatch count which will still allow a full match to be performed), palindrome clip threshold (how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment), and simple clip threshold (how accurate the match between any adapter sequence must be against a read); `LEADING:3` and `TRAILING:3` specifies minimum base quality at the ends of reads; `MINLEN:36` specifies the minimum read length for the read to be kept

Once we've filtered our reads, we can run Fastqc again if we'd like to check whether our data seem to have improved. Once we're satisfied, we can move on to aligning our reads to a genome.
