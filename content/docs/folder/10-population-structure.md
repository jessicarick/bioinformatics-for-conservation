---
title: 10 - Inferring Population Structure
type: docs
prev: docs/folder/9-estimating-ne
next: docs/folder/11-phylogenetics
---

## Inferring population structure or genetic clusters
Often, one of the first things that we need to assess about our data is whether the samples that we have sequenced comprise one panmictic population, or whether there is population structure present that we need to account for in our analyses. We may be interested in learning about population structure itself, or we may simply need to make sure that we're not violating an assumption of our analytical models by having samples that do not come from one cohesive population. Either way, we'll want to do similar analyses (although in the former case, we likely will want to do additional analyses based off of our findings).

### Visualizing structure: principal component analysis
The easiest, first-pass way of looking for population structure within our data usually is to conduct a principal component analysis (PCA) with all of our samples, which allows us to detect the greatest axes of divergence within our data and to visualize genetic relatedness among individuals. If we plot our our PCA results and see separate clusters, this will tell us that the population is *not* panmictic and that there are groups that are not randomly mating among themselves for some reason (and then we get to figure out what that reason is!).

I've already shared code for conducting a simple PCA on a `genlight` object using `dartR` in R, which can be found on the [Estimating HWE and PCAs](https://jessicarick.github.io/bioinformatics-for-conservation/docs/folder/8-hwe-pca/#plotting-a-pca-to-visualize-our-data) page. For the analyses below, you'll need the `AC1_variants_filtered.recode.vcf` file and related metadata (`AC1_ind_pops.txt`), which both can be found in our shared data folder on the HPC in the `week9_data/` directory.

```r
library(tidyverse)
library(dartR)

##################
## FIRST, PCA ####
##################
# import VCF
char_gl <- gl.read.vcf("AC1_variants_filtered.recode.vcf")

# perform pca (might take a bit because it's a big file)
char_gl_pca <- gl.pcoa(char_gl)

# plot pca
gl.pcoa.plot(char_gl_pca, char_gl)
gl.pcoa.plot(char_gl_pca, char_gl,
             xaxis=2, yaxis=3)
gl.pcoa.plot(char_gl_pca, char_gl,
             xaxis=3, yaxis=4)

# do these clusters correspond to sampling locations?
samplelist <- read_delim("data/AC1_ind_pops.txt", 
                         col_names=c("ind","pop"),
                         delim=" ")

char_gl@pop <- tibble(ind = indNames(char_gl)) %>%
  left_join(samplelist) %>%
  pull(pop) %>%
  as_factor()
  
gl.pcoa.plot(char_gl_pca, char_gl)

## yes, these clusters seem to be our sampling locations
## if not, then we could use k-means clustering to identify
## the groups first and assign inds to them, before performing DAPC

### CALCULATING FST
### How divergent are these groups?
char_fst <- reich.fst(char_gl, bootstrap=10, plot=TRUE, verbose=TRUE)

### RUNNING DAPC
### What if we're interested in the SNPs that contribute
### to separation among groups -- then DAPC
# (this also takes a while to run!)
char_dapc <- dapc(char_gl, n.pca=20, n.da=3)

scatter(char_dapc, scree.da = TRUE,
        bg="white", pch=20, cell=0,
        cstar=1, solid=0.4, cex=2, clab=0,
        leg=FALSE)

# loci on first axis
loadings1 <- char_dapc$var.contr[,1]
loadingplot(loadings1,
            cex.lab=0.5, srt=90, byfac=FALSE,
            xlab="Locus", threshold=quantile(loadings1, 0.995))

# loci on second axis
loadings2 <- char_dapc$var.contr[,2]
loadingplot(loadings2,
            cex.lab=0.5, srt=90, byfac=FALSE,
            xlab="Locus", threshold=quantile(loadings2, 0.995))
```

### Clustering analyses
To go beyond visualizing our data, we can use clustering algorithms to tell us which individuals, when grouped together, are in Hardy-Weinberg Equilibrium within their group but not among the groups. Two commonly used programs are [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) (slow for SNP-scale datasets, but has the related [FASTstructure](https://rajanil.github.io/fastStructure/) which is designed to better deal with them) and [ADMIXTURE](https://dalexander.github.io/admixture/) (similar to STRUCTURE, but faster and designed for SNP data), and other relatively well-used ones include [NGSAdmix](https://www.popgen.dk/software/index.php/NgsAdmix) (allows the use of genotype likelihoods; within the ANGSD ecosystem of programs), [entropy](https://bitbucket.org/buerklelab/mixedploidy-entropy/src/master/) (requires the use of genotype likelihoods, and can handle mixed ploidy data!), [fineSTRUCTURE](http://paintmychromosomes.com/) (which is designed for dense/high-resolution sequencing data), etc.

Together in class, we'll be running `ADMIXTURE` as an example, as it is relatively quick to run. We need our input data to be in PLINK format; I've already converted our `vcf` to the necessary PLINK files (`.bed`, `.bim`, and `.fam`). ADMIXTURE assumes that our SNPs are unlinked, so these SNPs have already been filtered to prune those that are in linkage disequilibrium with one another. We then need to choose the number of groups that we want to test in each run-- for these data, we'll run the program for K=1 to K=5. ADMIXTURE is installed as a module, so we can load it easily on the HPC.

```sh
module load admixture

# make directory for results and move into that directory
mkdir admix_results/
cd admix_results/

# run ADMIXTURE for k=1
admixture --cv /xdisk/jrick/consbioinf/shared_data/week9_data/AC1_variants_filtered.bed 1
```
This will give us results just for K=1. If we want to run the program for a bunch of K values, we can do so using a for loop:

```sh
# run ADMIXTURE for K=2 to K=5
for k in {2..5}; do
  echo "running for K=$k"
  admixture --cv /xdisk/jrick/consbioinf/shared_data/week9_data/AC1_variants_filtered.bed $k
done
```
Once we have all of our results, we will want to download all of the `.Q` files to our local computer so that we can visualize them in R. You'll also need to download the `admix_functions.R` script from the `week9_data/` directory on the HPC, which will give some custom functions for plotting ADMIXTURE results.

```r
library(tidyverse)
library(dartR)
source("scripts/admix_functions.R")

##########################################
## SECOND, PLOTTING ADMIXTURE RESULTS ####
##########################################

# load list of individuals and sampling sites
samplelist <- read_delim("data/AC1_ind_pops.txt", 
                         col_names=c("ind","pop"),
                         delim=" ")

# get list of output files from admixture
qfiles <- list.files(path = "data/", pattern = ".Q", full.names=TRUE)

# import all of the data
admix_results <- qfiles %>%
  map_dfr(., ~read_qfiles(.x)) 

# plot only k=2
plot_admix(admix_results)

# plot all of the k values
plot_admix(admix_results, min_k = 1, max_k=5)
```
