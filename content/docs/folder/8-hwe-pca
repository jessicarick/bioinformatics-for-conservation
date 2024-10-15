---
title: 8 - Estimating HWE and PCAs
type: docs
prev: docs/folder/7-relatedness-kinship-parentage
next: docs/folder/8-hwe-pca
---

This page is to provide some additional code that we didn't get to in class, which will be beneficial for completing Project 1.

## Estimating Hardy-Weinberg Equilibrium

There are two different methods for inferring whether sites are in Hardy-Weinberg equilibrium from our genomic data. One is using `VCFtools`, while the other ues the `dartR` package in R.

### Using VCFtools

To infer HWE for each site using VCFtools with your filtered VCF, you can use the following code:

```sh
vcftools --vcf variants_filtered.vcf --hardy --out vcf_stats/variants_filtered
```

This will output a p-value for each site from a Hardy-Weinberg Equilibrium test, as well as the observed numbers of homozygotes and heterozygotes and the corresponding expected numbers under HWE. These data can then be pulled into R for interpretation, using code such as that below:

```r
library(tidyverse)

## load and plot a histogram of the hwe file
hwe <- read_table("variants_filtered.hwe")

hwe %>%
  ggplot() +
  geom_histogram(aes(x=P_VALUE))

## summary of the hwe statistics
summary(hwe)
```

### Using dartR

To infer HWE for each site in R, we can use the `dartR` package with our filtered VCF file. If you have already loaded your VCF for performing other analyses, then you will not need to re-load it specifically for this step.

```r
library(tidyverse)
library(dartR)

char_gl <- dartR::gl.read.vcf("data/AC1_chr1_variants_filtered.recode.vcf")
char_gl

# calculate hwe per locus
# sig_only = TRUE means only report those that deviate from HWE
# sig_only = FALSE means include estimates for all snps
hwe <- gl.report.hwe(char_gl, sig_only = FALSE)
summary(hwe)

# plot histogram of hwe p-values
hwe %>%
  ggplot() +
  geom_histogram(aes(x=as.numeric(Prob)))

# plot heterozyosity vs hwe p-value
char_hwe %>%
  ggplot() +
  geom_point(aes(x=as.numeric(Het), y=as.numeric(Prob)))

```

## Plotting a PCA to visualize our data

[Principal component analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) is often used as a first-pass method for visualizing relationships among individuals in a genetic data set. For question 6 on your Project 1, you are asked to perform a PCA on your data. Although we haven't yet covered this in class, the following code will help you to perform and plot a quick PCA from your genlight object, again using `dartR`. Note that the calculation and plotting of the PCA may take a bit of time!

```r
library(tidyverse)
library(dartR)

char_gl <- dartR::gl.read.vcf("data/AC1_chr1_variants_filtered.recode.vcf")
char_gl

# perform pca
char_gl_pca <- gl.pcoa(char_gl)

# plot pca
gl.pcoa.plot(char_gl_pca, char_gl)

# alternative plot, if we want to plot it ourselves
char_gl_pca$scores %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point() +
  theme_bw() +
  xlab(paste0("PC1, ", round(char_gl_pca$eig[1]/sum(char_gl_pca$eig)*100,1), "% of variation"))  +
  ylab(paste0("PC2, ", round(char_gl_pca$eig[2]/sum(char_gl_pca$eig)*100,1), "% of variation"))

```

In our PCA plots, each point represents one individual in our dataset, and individuals that are closer to one another on the plot are more genetically similar to one another. If we see clusters of points, then this suggests that we have clusters of individuals in our dataset that have reduced gene flow between those clusters. 
