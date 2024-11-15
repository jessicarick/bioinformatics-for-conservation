---
title: 12 - Landscape genomics
type: docs
prev: docs/folder/11-phylogenetics
next: docs/folder/13-hybridization
---

## Landscape genetics and landscape genomics
Generally, landscape genetics/genomics combines environmental data with genetic data to learn something about the organism that we're interested in. This kind of research combines tools from population genetics and landscape ecology, and has two different general flavors:

(1) **Connectivity analyses** (generally fall under "landscape genetics"), where we're interested in the influence of environmental features on the genetic structure or gene flow between populations. In these analyses, the goal is to identify landscape features that enhance or restrict gene flow using what we call resistance surfaces. In these, each landscape variable is assigned a resistivity (i.e., how resistant we think it is to gene flow in the given population) and then these resistance surfaces are used to calculate effective distances between individuals (or populations). These "landscape distances" are then compared to genetic distances in a formal modeling framework to determine the best fit to our genetic data (high correlation between landscape distance and genetic distance for a given landscape variable or set of variables = better fit).

(2) **Genotype-environment associations** or outlier tests (generally fall under "landscape genomics"), where we're interested in identifying the environmental factors that have shaped present-day adaptive variation and the genetic variants that drive local adaptation. These are the type of analyses that the [Rellstab et al. (2015)](https://doi.org/0.1111/mec.13322) paper focuses on, and that we'll focus on as well for the coding exercises in this class. In these analyses, we're interested in looking for specific loci under natural selection in specific environments, and associating those with specific environmental variables.

## Genome scans and outlier analyses
There are many different flavors of genome scans/outlier analyses/genotype-environment association tests, depending on the type of data that you have and the question that you're interested in. Generally, they involve testing each locus or SNP in your dataset for a correlation between allele frequencies at that site (or in that window) and one or more environmental variables of interest, while accounting for background population structure. These types of analyses can have high rates of false positives, and the results should be interpreted with appropriate caution! The most robust method is to use mutiple different types of analyses (with different underlying assumptions), and then to look for overlap among the loci identified.

For our in-class purposes, we are going to be performing a [latent factor mixed model](https://doi.org/10.1093/molbev/mst063) analysis on a publicly-available dataset of sunflower (*Heliantus annus*) genotypes that were generated in XX (NCBI SRA Project XX). To perform the LFMM analysis, we'll be using the R package `LEA`. The data need to be in the `012` format, where 0=homozygous for the dominant allele, 1=heterozygous, and 2=homozygous for the recessive allele. We can convert a normal VCF into this format using VCFtools (e.g., `vcftools --vcf variants_filtered.vcf --012 --out variants_filtered`). Missing data also need to be coded as a `9`. I've already converted the genomic data into the appropriate format, and the ready-to-go data (named `hannus.freebayes.80.1.lfmm`) can be found in the `week11_data/` directory within our `shared_data/` folder on the UA HPC. There are also associated metadata available in the `hannus_selected_inds_metadat.csv` file within the same directory.

In addition to our genomic data, we'll also need data for whatever environmental/landscape variable we're interested in. For today, we'll be looking for SNP associations with precipitation, and so we'll be using the [WorldClim](https://www.worldclim.org/data/index.html) dataset for Annual Precipitation (part of the [BioClim](https://www.worldclim.org/data/bioclim.html) dataset). These data are also available to you in the `week11_data/` directory, in the file named `wc2.1_10m_bio_12.tif`.

In class, we used the following R code (note that I will add more explanations of this code later). Note that this code has been updated to no longer use the deprecated `rgdal` and `raster` packages, but instead now uses `terra` and `sp`!

```r
# lfmm analyses using LEA
# week 11 -- landscape genomics

#BiocManager::install("LEA")
#BiocManager::install("qvalue")
#install.packages("sp") 
#install.packages("terra")

library(tidyverse)
library(LEA)
library(qvalue)
library(sp)
library(terra)

### GENETIC DATA
# import .lfmm file
dat_lfmm <- read_table("data/hannus.freebayes.80.1.lfmm",
                       col_names=FALSE)
write.lfmm(dat_lfmm,"data/hannus.lfmm")

## Infer individual admixture coefficients using snmf
# (to be used to account for neutral variation)
# main options 
# K = number of ancestral populations 
# entropy = TRUE computes the cross-entropy criterion, 
geno_snmf <- NULL
geno_snmf <- snmf("data/hannus.lfmm", 
               K = 1:6, 
               entropy = TRUE, 
               repetitions = 10,
               project = "force",
               CPU=1)
               
# plot cross-entropy criterion for all runs in the snmf project 
# minimum = best fit
plot(geno_snmf, col = "blue", pch = 19, cex = 1.2)

# here, it looks like K=2 is the best fit for the data

## CLIMATE DATA
#Load sample names and coordinates
sample.coord <- read_csv("data/hannus_selected_inds_metadat.csv")
sample.coord <- sample.coord %>%
  mutate(Latitude = as.numeric(Lat),
         Longitude = as.numeric(Long))

#Define the spatial projection system that the points are in (usually WGS84)
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  
sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('Longitude','Latitude')], 
                                          proj4string=CRS(crs.wgs), data=sample.coord)
sample.coord.sv <- vect(sample.coord.sp)

# import climate data (from WorldClim)
# here, we're using the bioClim variable 12 (annual precipitation)
clim.layer <- rast("data/wc2.1_10m_bio_12.tif")

#Extract the climate data for each point (projection of climate layer and coordinates must match)
clim.points <- extract(clim.layer, sample.coord.sv) 

#Combine the sample coordinates with the climate data points and save for use with GF tomorrow
clim.points <- cbind(sample.coord, clim.points)  
write.table(clim.points, "data/clim.points", 
            sep="\t", quote=F, row.names=F)  
clim.points 

#Save climate data without sample names, latitude, longitude, or column names for LFMM
clim.env <- clim.points$wc2.1_10m_bio_12
colnames(clim.env) <- NULL
clim.env
write.table(clim.env, "data/clim.env", 
            sep="\t", quote=F, row.names=F, col.names=F) 

## Now, running LFMM
# main options: 
# K = the number of latent factors
# Runs with K = 2 using 3 repetitions.
# in a real run, you'll want to run longer than 5000 iterations!
geno_lfmm <- NULL
geno_lfmm <- lfmm("data/hannus.lfmm", 
               "data/clim.env", 
               K = 2, 
               repetitions = 3, 
               CPU = 1,
               iterations = 5000,
               burnin = 500, 
               project = "new")

# now, we can combine our three runs to get new p-values
# by taking the median of the z scores
geno_lfmm_z <- z.scores(geno_lfmm, K = 2, d = 1) # here, K needs to match your chosen K above
geno_lfmm_z <- apply(geno_lfmm_z, 1, median)

# next, we have to use z-scores to calculate lambda 
# ("genomic inflation factor")
geno_lfmm_lambda <- median(geno_lfmm_z^2)/qchisq(0.5, df = 1)
geno_lfmm_lambda

# use lambda to calculate adjusted p-values
geno_lfmm_p_adj <- pchisq(geno_lfmm_z^2/geno_lfmm_lambda, 
                         df = 1, lower = FALSE)

# finally, correct p-values for multiple testing
geno_lfmm_q_final <- qvalue(geno_lfmm_p_adj)$qvalues

# how many are outliers at p < 0.05?
sum(geno_lfmm_q_final < 0.05)

# manhattan plot of q values!
plot(-log10(geno_lfmm_q_final), 
     pch = 19, 
     col = "blue", 
     cex = 0.5, 
     xlab = "SNP")

# add lines for p < 0.05 and p < 0.005 significance thresholds
abline(h=-log10(0.05))
abline(h=-log10(0.005), lty=3)
```
This should give you a Manhattan plot that looks like the following, where each dot represents one SNP in our dataset, ordered by their position in the genome, and dots with larger y-values are ones that have a stronger correlation with our precipitation data! The solid line is our p < 0.05 significance threshold and the dotted line shows where p < 0.005.

![Figure showing a Manhattan plot of LFMM results.](/content/img/lfmm_manhattan_example.png)
