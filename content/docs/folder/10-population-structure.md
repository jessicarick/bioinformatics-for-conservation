---
title: 10 - Inferring Population Structure
type: docs
prev: docs/folder/9-estimating-ne
next: docs/folder/10-population-structure
---

## Inferring population structure or genetic clusters
**Please note that this page is a work in progress

Often, one of the first things that we need to assess about our data is whether the samples that we have sequenced comprise one panmictic population, or whether there is population structure present that we need to account for in our analyses. We may be interested in learning about population structure itself, or we may simply need to make sure that we're not violating an assumption of our analytical models by having samples that do not come from one cohesive population. Either way, we'll want to do similar analyses (although in the former case, we likely will want to do additional analyses based off of our findings).

### Visualizing structure: principal component analysis
The easiest, first-pass way of looking for population structure within our data usually is to conduct a principal component analysis (PCA) with all of our samples, which allows us to detect the greatest axes of divergence within our data and to visualize genetic relatedness among individuals. If we plot our our PCA results and see separate clusters, this will tell us that the population is *not* panmictic and that there are groups that are not randomly mating among themselves for some reason (and then we get to figure out what that reason is!).

I've already shared code for conducting a simple PCA on a `genlight` object using `dartR` in R, which can be found on the [Estimating HWE and PCAs](https://jessicarick.github.io/bioinformatics-for-conservation/docs/folder/8-hwe-pca/#plotting-a-pca-to-visualize-our-data) page.

### Clustering analyses
To go beyond visualizing our data, we can use clustering algorithms to tell us which individuals, when grouped together, are in Hardy-Weinberg Equilibrium within their group but not among the groups. Two commonly used programs are [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) (slow for SNP-scale datasets, but has the related [FASTstructure](https://rajanil.github.io/fastStructure/) which is designed to better deal with them) and [ADMIXTURE](https://dalexander.github.io/admixture/) (similar to STRUCTURE, but faster and designed for SNP data), and other relatively well-used ones include [NGSAdmix](https://www.popgen.dk/software/index.php/NgsAdmix) (allows the use of genotype likelihoods; within the ANGSD ecosystem of programs), [entropy](https://bitbucket.org/buerklelab/mixedploidy-entropy/src/master/) (requires the use of genotype likelihoods, and can handle mixed ploidy data!), [fineSTRUCTURE](http://paintmychromosomes.com/) (which is designed for dense/high-resolution sequencing data), etc.
