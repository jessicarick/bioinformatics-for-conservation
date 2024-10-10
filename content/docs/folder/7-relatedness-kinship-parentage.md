---
title: 7 - Estimating Relatedness, Kinship, and Parentage
type: docs
prev: docs/folder/6-estimating-diversity
next: docs/folder/7-relatedness-kinship-parentage
---

## Estimating relatedness, kinship, and parentage
All three of these tasks are somewhat related, but can be approached in several different ways. Essentially, they all entail quantifying the relationships between individuals in a dataset, and then in the case of kinship and parentage, then binning those relationships into parent-offpring, full sibling, half sibling, etc. according to theory on how related individuals should be for each of these categories. The most recent methods for these sorts of estimation methods use Bayesian or maximum likelihood inference to quantify our uncertainty in a given estimate of kinship or a given relationship between individuals. 

### Relatedness and kinship using ngsRelate
For our class, we're going to approach this two different ways. First, we're going to estimate relatedness among all individuals in a dataset using [NgsRelate](https://github.com/ANGSD/NgsRelate), which is part of the ANGSD ecosystem of programs. NgsRelate will actually co-estimate both relateness and inbreeding coefficients, so we can efficiently get estimates for both of these measures using the same program. We can start analyses in ngsRelate from either `.bam` or `.vcf` files (or from a `.glf` file output from ANGSD), depending on our analysis pipeline. When using a `.vcf` file, ngsRelate will make calculations based on genotype likelihoods, but we can also specify that we would prefer it to use genotype calls if we don't have likelihoods.

For the simplest estimation of relatedness among individuals, we can use the following code:

```sh
# save path to program installation
alias ngsrelate="/xdisk/jrick/programs/ngsRelate/ngsRelate

# note that your VCF file needs to be zipped here! we'll first do that using bgzip
bgzip -c AC1_chr1_variants_filtered.recode.vcf > AC1_chr1_variants_filtered.vcf.gz

# run the program using VCF as input
# also feed it a list of individual IDs
ngsRelate -h AC1_chr1_variants_filtered.vcf.gz -z ind_list_AC1.txt -O AC1_chr1_variants_filtered.rel
```

From the output, you'll notice that ngsRelate has a default minor allele frequency filter of 0.05. We'll leave this be for now, but could adjust that number in the future. The output from ngsRelate has individual IDs in the first two columns, and then for each pair of individuals has a number of different calculated values. You can read about the output format on the ngsRelate GitHub page here: [https://github.com/ANGSD/NgsRelate?tab=readme-ov-file#output-format](https://github.com/ANGSD/NgsRelate?tab=readme-ov-file#output-format). Two of the measures that we may be interested in are called `rab` (column 13; pairwise relatedness from [Hedrick et al. 2015](https://academic.oup.com/jhered/article/106/1/20/2961876)) and `KING` (column 31; kinship from [Waples et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/mec.14954)). To make sense of these estimates, we'll probably want to pull them into R and plot their distributions, and then we can look closer at specific individuals that may be more related than expected. Generally, kinship values of 0.5 indicate first-degree relatives (full siblings or parent-offspring relationships), while a kinship value of 1.0 indicates self or monozygotic twins.

```r
library(tidyverse)

relatedness <- read_table("AC1_chr1_variants_filtered.rel")

# distribution of kinship coefficients
relatedness %>%
  ggplot(aes(x=KING)) +
  geom_histogram()

# relationship between rab relatedness and KING coefficient
relatedness %>%
  ggplot(aes(x=rab,y=KING)) +
  geom_point()
```

### Parentage analysis 
Parentage analysis is a little bit trickier than just estimating relatedness or kinship among individuals. However, it can be a really important step in some conservation-related genomic analyses. Again, there are several programs and R packages that can be used to infer parentage; some common ones include [the Sequoia R package](https://jiscah.github.io/articles/vignette_main/book/index.html), [COLONY](https://www.zsl.org/about-zsl/resources/software/colony), [CERVUS](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp) (more commonly used for microsatellite data), and the [CKMRsim R package](https://eriqande.github.io/CKMRsim/).
