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
