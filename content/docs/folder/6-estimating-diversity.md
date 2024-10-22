---
title: 06 - Estimating Genetic Diversity
type: docs
prev: docs/folder/5-variant-filtering
next: docs/folder/7-relatedness-kinship-parentage
---

## Estimating genetic diversity
There are several different measures that are all considered estimates of "genetic diversity" in a population, and numerous different programs that can be used for calculating each of these. For the purposes of our class, we're going to focus on the following types of genetic diversity:
* **Observed heterozygosity**: a measure of the mean rate of per-individual heterozygosity in the population; usually calculated by dividing the number of heterozygous SNPs by the total number of SNPs for each individual, and then taking the average of this; can be biased by sample size
* **Expected heterozygosity** (H<sub>E</sub> or H<sub>S</sub>): a measure of the probability that a pair of randomly sampled allele copies from a population will be different from one another; calculated as 1 - the sum of the squared allele frequencies at each locus in each individual (comes from the Hardy-Weinberg equation $p^2 + 2pq + p^2 = 1$)
* **Watterson's theta** ($\theta~W$): an estimate of the number of segretating (polymorphic) sites within a group; related to $N_e$ via the equation $\theta = 4N_e\mu$, where $\mu$ is the mutation rate
* **Nucleotide diversity** ($\pi$): an estimate of the average number of differences (mismatches) between any pair of DNA sequences in a group
* **(Tajima's D)**: not an estimate of diversity *per se*, but rather the ratio of $\theta~W$ and $\pi$ in a population; under neutral expectations for a population at equilibrium, we expect the two to be equal to one another, and deviations from this expectation provide evidence for an excess of rare alleles (when negative) or a deficit of rare alleles (when positive) compared to expectations
* **Runs of homozygosity**: a measure of the length and frequency of regions of the genome that are autozygous (i.e., identical by descent from a common ancestor), whcih can be used to quantify levels of inbreeding and the genetic relatedness between individuals

To calculate these measures of diversity, we'll be using two different methods: the program [ANGSD](https://www.popgen.dk/angsd/index.php/ANGSD#Overview), which requires sorted `bam` files as input, and R, where we can work directly from our filtere `VCF` that we created.

## Genetic diversity using ANGSD
The following code can be used to generate estimates of genetic diversity from `bam` files using `ANGSD`. This program is not available as a module on the UA HPC, so we first have to specify where it is installed so that the computer knows where to find it. The method in ANGSD first requires estimating the site allele frequency likelihood (SAF; an estimate of the allele frequencies within the population), followed by an estimation of the site frequency spectrum (SFS), which is then used for estimating overall genetic diversity.

```sh
# these programs aren't available as modules,
# so we need to indicate where to find them
alias angsd="/xdisk/jrick/programs/angsd/angsd"
alias realSFS="/xdisk/jrick/programs/angsd/misc/realSFS"
alias thetaStat="/xdisk/jrick/programs/angsd/misc/thetaStat"

# create an output directory
mkdir angsd_out

# create SAF file for chromosome 1 with 10 individuals (runs relatively quickly)
# all comments will need to be removed and flags will all have to be
# on the same line for code to work
angsd \
	-b /xdisk/jrick/consbioinf/shared_data/char_bamlist_short.txt \ # list of bamfiles to include
	-ref /xdisk/jrick/consbioinf/shared_data/char_reference/Salvelinus_spp_genome.fasta  \ # reference genome
	-anc /xdisk/jrick/consbioinf/shared_data/char_reference/Salvelinus_spp_genome.fasta  \ # reference genome
	-out angsd_out/char_bamlist_short \ # output name
	-uniqueOnly 1 \ # keep only reads with 1 best mapping location
	-minMapQ 20 \ # minimum mapping quality of 20
	-minInd 2 \ # sites have to have at least 2 individuals with genotype calls
	-setMinDepth 10 \ # keep only if total read depth across individuals is > 10
	-setMaxDepth 1000 \ # keep only if total read depth across individuals is < 1000
	-doCounts 1 
	-GL 1 \ # use samtools model for GLs 
	-doSaf 1 \ # Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
	-r NC_036858.1  # only working with chromosome 1

# then estimate the SFS from the SAF
realSFS char_bamlist_short.saf.idx > char_bamlist_short.sfs

# now, we use the sfs to estimate diversity statistics for the population
realSFS saf2theta char_bamlist_short.saf.idx -sfs char_bamlist_short.sfs -outname char_bamlist_short

# and finally calculate the genome-wide summary statistics
thetaStat do_stat char_bamlist_short.thetas.idx 

# in the output file, tW = watterson's theta; tP = pi; Tajima = Tajima's D
# the first two will need to be divided by the number of sites to get proportional estimates
# these can also be estimated in windows if that is desired

```

## Genetic diversity using R
To generate estimates using `R`, we can use the following code, which depends heavily on the [dartR package](https://green-striped-gecko.github.io/dartR/). If you do not already have `dartR` installed, then you can install it using `install.packages("dartR")` followed by `gl.install.vanilla.dartR()`. If you run into errors during the installation, there is some helpful additional installation information found at [https://github.com/green-striped-gecko/dartR/wiki/Installation-tutorial#installation](https://github.com/green-striped-gecko/dartR/wiki/Installation-tutorial#installation).

```r
library(dartR) 
library(adegenet) 
library(tidyverse)

## import vcf and convert to genlight object
char_gl <- dartR::gl.read.vcf("data/AC1_chr1_variants_filtered.recode.vcf")
char_gl

# how many individuals do we have?
nInd(char_gl)

# how many snps?
nLoc(char_gl)

# calculate missing data per locus -- do we need to do more filtering?
gl.report.callrate(char_gl,method="loc")

# calculate missing data per individual -- do we need to remove anyone?
char_miss_ind <- rowSums(is.na(as.matrix(char_gl)))/nLoc(char_gl)

# calculate heterozygosity per individual
char_het <- gl.report.heterozygosity(char_gl, method="ind")

# look at missing data vs heterozygosity
char_het %>%
  mutate(missing=char_miss_ind) %>%
  ggplot(aes(x=missing, y=Ho)) +
  geom_point()

# perhaps we should filter individuals for >80% missing data as well
char_gl_filt <- gl.filter.callrate(char_gl, method="ind",
                                   threshold=0.2, mono.rm=TRUE,
                                   recalc=TRUE)

# now, recalculate missing data and heterozygosity
# calculate missing data per individual -- do we need to remove anyone?
char_miss_ind_filt <- rowSums(is.na(as.matrix(char_gl_filt)))/nLoc(char_gl_filt)
char_het_filt <- gl.report.heterozygosity(char_gl_filt, method="ind")

char_het_filt %>%
  mutate(missing=char_miss_ind_filt) %>%
  ggplot(aes(x=missing, y=Ho)) +
  geom_point()

# this looks a little cleaner now, so let's calculate diversity metrics!
char_div <- gl.basic.stats(char_gl_filt)
char_div
```

