---
title: 5 - Variant filtering
type: docs
prev: docs/folder/4-variant-calling
next: docs/folder/5-variant-filtering
---

## Variant Filtering
Now that we have our raw variants file, we are going to focus on how to filter these variants. While this may seem like a trivial task, it is exceedingly important and can influence the results of your downstream analyses (add links to papers), so should be thought about carefully.

The first thing that we may want to check is to see how many variants we have in the unfiltered `.vcf` file. To do this, we can use `bcftools` again, with the `view` command:

```sh
bcftools view -H variants.vcf | wc -l
```
We will likely have many variants in this unfiltered file, which is great! However, many of them are likely not useful for analyses because they only have data for a few individuals, they don't have sufficient read depth, or they are low quality. To get a better idea of how we might want to filter our variants, we can generate statistics on our `.vcf` file using `vcftools`. In particular, we are going to generate some statistics about depth, quality, minor allele frequency, and missing data.

```sh
# report mean depth per individual
vcftools --vcf variants.vcf --depth --out vcf_stats/variants

# report mean depth per site
vcftools --vcf variants.vcf --site-mean-depth --out vcf_stats/variants

# report per-individual missingness
vcftools --vcf variants.vcf --missing-indv --out vcf_stats/variants

# report per-site missingness
vcftools --vcf variants.vcf --missing-site --out vcf_stats/variants

# report per-individual heterozyosity and Fis
vcftools --vcf variants.vcf --het --out vcf_stats/variants
```
Once we have each of these reports, we can use `R` to summarize and visualize their distributions to inform our further filtering of the data. To do this, we can download the output files (all contained within `vcf_stats/`) using `scp` on the command line or the GUI interface at [https://ood.hpc.arizona.edu](https://ood.hpc.arizona.edu) and use a script such as the following code:

```r
## script for looking at vcftools stats output
## WFSC 496/596 in-class week 5

## load necessary packages
library(tidyverse)

prefix="AC1_chr1_variants_bcftools"

## load and plot idepth file
idepth <- read_table(paste0(prefix,'.idepth'))

idepth %>%
  ggplot() +
  geom_histogram(aes(x=MEAN_DEPTH))

## load and plot imiss file
imiss <- read_table(paste0(prefix,'.imiss'))

imiss %>%
  ggplot() +
  geom_histogram(aes(x=F_MISS))

## load and plot het file
ihet <- read_table(paste0(prefix,'.het'))

ihet %>%
  ggplot(aes(x=`O(HOM)`)) +
  geom_histogram()

# transform into observed heterozygosity
ihet %>%
  mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
  ggplot(aes(x=HET_O)) +
  geom_histogram()

# what about inbreeding (F)?
ihet %>%
  ggplot(aes(x=F)) +
  geom_histogram()

## we can combine to look at depth vs missing data
idepth %>%
  left_join(imiss,by="INDV") %>%
  ggplot(aes(x=MEAN_DEPTH,y=F_MISS)) +
  geom_point()

# who is that weird outlier?
ihet %>%
  mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
  filter(HET_O > 0.3)

## can also combine to look at het vs missing data
imiss %>%
  left_join(ihet,by="INDV") %>%
  mutate(HET_O = 1-(`O(HOM)`/N_SITES)) %>%
  ggplot(aes(x=F_MISS,y=HET_O)) +
  geom_point()

## load ldepth file
ldepth <- read_table(paste0(prefix,'.ldepth.mean'))

ldepth %>% 
  ggplot() +
  geom_histogram(aes(x=MEAN_DEPTH))

# histogram not super informative, maybe a plot across position will be better?
ldepth %>% 
  ggplot() +
  geom_point(aes(x=POS,y=MEAN_DEPTH))

summary(ldepth$MEAN_DEPTH)

## load lmiss file
lmiss <- read_table(paste0(prefix,'.lmiss'))

lmiss %>% 
  ggplot() +
  geom_histogram(aes(x=F_MISS))

# maybe a plot across position will be also informative?
lmiss %>% 
  ggplot() +
  geom_point(aes(x=POS,y=F_MISS))

summary(lmiss$F_MISS)

```

