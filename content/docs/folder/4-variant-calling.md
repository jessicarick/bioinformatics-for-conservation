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

If you are interested in reading more about variant calling and the different options, [Nielsen et al. (2011)](https://doi.org/10.1038/nrg2986) provides a nice review, with a focus on mitigating uncertainty in genotype calling.

### Brief aside: genotype likelihoods and genotype probabilites
Often, we are not 100% certain about a base call at a specific site, especially when we have lower read depth. In these situations, we can use genotype probabilities or likelihoods instead of genotype calls, which allow us to quantify our uncertainty (or certainty) in a given genotype and to use that in downstream calculations. For moderate or low sequencing depths, genotype calling based on fixed cutoffs will typically lead to under-calling of heterozygous genotypes. The advantages of the probabilistic methods are that they provide measures of statistical uncertainty when calling genotypes, they lead to higher accuracy of genotype calling, and they provide a natural framework for incorporating information regarding allele frequencies and patterns of LD into our inferences.

Genotype likelihoods are the probability of the observed data (i.e., the reads) given a specific genotype. The likelihoods are presented in threes for diploid data, which corresponds to the likelihood that the true genotype is AA, AB, and BB. All four of the variant callers above can provide genotype likelihoods in the output VCF files. Many downstream methods require called genotypes rather than likelihoods; however, there are more and more methods being developed that will do analyses on likelihoods, so those are beneficial to know about. While likelihoods are necessary for low-depth data, it is also good practice to use them where possible even for higher-depth data, as they allow you to incorporate uncertainty into your downstream analyses instead of ignoring this potential source of error. There is a more in-depth explanation of genotype likelihoods for those interested in the variant calling section of [Eric Anderson's Bioinformatics Handbook](https://eriqande.github.io/eca-bioinf-handbook/variant-calling.html).

## Calling variants using BCFtools and SAMtools
For today, we'll be working with the BCFtools/SAMtools pipeline for calling variants and dealing with "hard" genotype calls, as detailed [here](https://samtools.github.io/bcftools/howtos/variant-calling.html). To do this, we'll want to have all of the sorted `.bam` files in the same directory, and we'll also need to know the path to our reference genome. Once we have both of those, we can use the following code to call variants across all of our samples combined:

```sh
module load bcftools
bcftools mpileup -Ou -f reference.fa bamfiles/*.sorted.bam | bcftools call -mv -Ou -o variants.vcf
```

To see what each of the flags that we use in this code means (and to see what other options exist), we can check out the BCFtools `mpileup` and `call` manual pages ([https://samtools.github.io/bcftools/bcftools.html#mpileup](https://samtools.github.io/bcftools/bcftools.html#mpileup) and [https://samtools.github.io/bcftools/bcftools.html#call](https://samtools.github.io/bcftools/bcftools.html#call)). For example, if we run this code, we will see a warning that only the first 250 reads for each region are examined and we may want to increase this number. To increase the read depth analyzed, we can add a `-d 1000` option to increase it to the first 1000 reads. We also likely want to retain read depth information (so that we can filter based on depth later on), in which case we may want to add the `-a AD,DP` flag after our `bcftools mpileup` command.
