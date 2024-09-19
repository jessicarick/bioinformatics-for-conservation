---
title: Quality Assessment and Filtering
type: docs
prev: docs/folder/prevariant-filtering
next: 
---

## Aligning reads to a genome
Once the reads look good, then we need to align them to a genome so that we can make sense of them. For this, we need (1) a program to use, (2) the reads we want to align, and (3) the genome that we want to align the reads to. The common alignment programs that are used for Illumina short reads include [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) and [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). These are both available as modules: `module load bwa` // `module load bowtie2`.

With your own data in the future, you may find that other aligners may be better for your purposes, depending on what kind of data you’re working with. There also are “workflow” programs that can be used to do many of these steps together: [STACKS](https://catchenlab.life.illinois.edu/stacks/manual/), [ipyrad](https://ipyrad.readthedocs.io/en/master/index.html), etc. (but these usually use bwa or bowtie2 under the hood). These can be nice because they don't require as much hands-on work, but they can also sometimes become "black boxes" that can get you into trouble if you're not thinking explicitly about the settings that you're using.

### Preparing the reference genome
We also need a reference genome to use! Look on NCBI genome database if you don’t have one, you can download the genome file (for our in-class work, we’ll be using the Arctic charr reference genome: [https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002910315.2/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002910315.2/), which I’ve already downloaded to our shared directory). I’ve already indexed the genome, but to do so with your own genome, you’ll need to use `bwa index genome.fa.gz` and `bowtie2-build genome.fa.gz genome_prefix` (these can take a while, depending on the size of your genome). Indices are necessary to allow the programs to more quickly access the genome once they start trying to align reads.


