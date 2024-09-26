---
title: 5 - Variant filtering
type: docs
prev: docs/folder/4-variant-calling
next: docs/folder/5-variant-filtering
---

** Note that this page is not yet complete **

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
Once we have each of these reports, we can use `R` to summarize and visualize their distributions to inform our further filtering of the data.

