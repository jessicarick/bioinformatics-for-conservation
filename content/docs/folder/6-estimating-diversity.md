---
title: 6 - Estimating Genetic Diversity
type: docs
prev: docs/folder/5-variant-filtering
next: docs/folder/6-estimating-diversity
---

## Estimating genetic diversity
There are several different measures that are all considered estimates of "genetic diversity" in a population, and numerous different programs that can be used for calculating each of these. For the purposes of our class, we're going to focus on the following types of genetic diversity:
* Observed heterozygosity: a measure of the mean rate of per-individual heterozygosity in the population; usually calculated by dividing the number of heterozygous SNPs by the total number of SNPs for each individual, and then taking the average of this; can be biased by sample size
* Expected heterozygosity (H<sub>E</sub> or H<sub>S</sub>): a measure of the probability that a pair of randomly sampled allele copies from a population will be different from one another; calculated as 1 - the sum of the squared allele frequencies at each locus in each individual (comes from the Hardy-Weinberg equation $p^2 + 2pq + p^2 = 1$)
* Watterson's theta ($\theta$<sub>w</sub>): an estimate of the number of segretating (polymorphic) sites within a group; related to $N<sub>e</sub>$ via the equation $\theta = 4N<sub>e</sub>\mu$, where $\mu$ is the mutation rate
* Nucleotide diversity ($\pi$)
* (Tajima's D)
* Runs of homozygosity
