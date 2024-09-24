---
title: 4 - Variant calling
type: docs
prev: docs/folder/3-read-alignment
next: docs/folder/4-variant-calling
---

** Note that this page is not yet complete **

## Variant Calling Software
Once we have our `.bam` files for all of our individuals, we can move forward with variant calling. There are many, many variant callers out there, and which one you want to use will depend on the specifics of your data and your planned analyses with those data. Some of the most popular variant callers for short read sequencing data and their features include:
* [GATK](https://gatk.broadinstitute.org/hc/en-us) -- a very common analysis pipeline developed by the Broad Institute (i.e., developed for use in humans) that includes the `HaplotypeCaller` module for calling variants within a population 
* [SAMtools/BCFtools](https://www.htslib.org/)
* [ANGSD](https://www.popgen.dk/angsd/index.php/ANGSD#Overview) -- a software for analyzing next generation sequencing data that has a lot of flexibility and can handle a number of different input types. Most of its methods take genotype uncertainty into account instead of basing the analysis on called genotypes, which is especially useful for low and medium depth data.
* [Freebayes](https://github.com/freebayes/freebayes) -- a Bayesian variant caller that uses a different method of variant detection (i.e., haplotype-based variant detection) than other popular software (see their documentation for more details)

## Calling variants using BCFtools and SAMtools
For today, we'll be working with the BCFtools/SAMtools pipeline for calling variants, as detailed [here](https://samtools.github.io/bcftools/howtos/variant-calling.html). To do this, we'll want to have all of the sorted `.bam` files in the same directory, and we'll also need to know the path to our reference genome. Once we have both of those, we can use the following code to call variants across all of our samples combined:

```sh
module load bcftools
bcftools mpileup -Ou -f reference.fa bamfiles/*.sorted.bam | bcftools call -mv -Ou -o variants.vcf
```
To see what each of the flags that we use in this code means (and to see what other options exist), we can check out the BCFtools `mpileup` and `call` manual pages ([https://samtools.github.io/bcftools/bcftools.html#mpileup](https://samtools.github.io/bcftools/bcftools.html#mpileup) and [https://samtools.github.io/bcftools/bcftools.html#call](https://samtools.github.io/bcftools/bcftools.html#call)). For example, if we run this code, we will see a warning that only the first 250 reads for each region are examined and we may want to increase this number. To increase the read depth analyzed, we can add a `-d 1000` option to increase it to the first 1000 reads. We also likely want to retain read depth information (so that we can filter based on depth later on), in which case we may want to add the `-a AD,DP` flag after our `bcftools mpileup` command.
