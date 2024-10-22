---
title: 09 - Estimating Effective Population Size
type: docs
prev: docs/folder/8-hwe-pca
next: docs/folder/10-population-structure
---

## Effective population size
The slides that contain the material covered in class about effective population size can be found [here](https://drive.google.com/file/d/17YW8fUv45wZc31xmsug2CBWeR_z8Hw6v/view?usp=sharing). The Zoom recording of our class covering this material can be watched or downloaded [here](https://drive.google.com/file/d/1UZ4M4UrU4fTzckfhJEZ_l8htZfmnh9gp/view?usp=sharing).

## Estimating effective population size
As with many of the different things that we're doing in class, there are many different methods available for estimating Ne, dependent on your data and the assumptions you're comfortable with making about your system. For our purposes, we're going to use the linkage disequilibirum method as implemented in [NeEstimator](), using `dartR` to run the program and process the results from within `R`. 

First, you'll need to download the `NeEstimator` program to your computer (from [https://www.molecularfisherieslaboratory.com.au/download-software/ ](https://github.com/user-attachments/assets/fe06427d-a837-4f2b-b789-c9f030d5434d), where you'll click on the first blue button), unzip it, and make note of where the folder is located on your computer. Once you have that, you'll also need a filtered VCF with variants that are only biallelic SNPs. In addition, your individual names and locus names within your VCF will need to not include any decimal places. We used the "clean" VCF produced from the char data on Chromosome 1, located at `/xdisk/jrick/consbioinf/shared_data/char_vcf/AC1_chr1_variants_clean_filtered.recode.vcf`, which you'll need to download to the `data/` directory within your R project. Once we do that, we can run the following code:

```r
library(tidyverse)
library(dartR)

# specify path to NeEstimator folder
ne_path <- "~/Downloads/Zip Folder_32_bit_191125" # this needs to match where your NeEstimator ended up

# import clean char VCF
char_gl <- gl.read.vcf("data/AC1_chr1_variants_clean_filtered.recode.vcf")
char_gl

# calculate Ne from VCF
char_ne <- gl.LDNe(char_gl,
                    outfile = "charLD.txt",
                    neest.path = ne_path,
                    critical = c(0, 0.05), # give results for using both MAF<0 and MAF<0.05
                    singleton.rm = TRUE, # remove any alleles only found in one individual
                    mating = "random") # mating system, other option is "monogamy"
char_ne
```
This will produce a plot with the Ne estimates (and confidence intervals) using three different filtering schemes for sites: (1) removing no sites (MAF = 0), (2) removing any sites where the minor allele frequency is less than 0.05 (MAF = 0.05), and (3) removing any sites where the minor allele is only found in one individual. In small datasets, including sites with low frequency alleles can bias your results (see discussion and math in [Waples & Do 2010](https://doi.org/10.1111/j.1752-4571.2009.00104.x)); however, this effect seems to be minimal in larger datasets and isn't generally a concern when working with SNP data. In the output of the `char_ne` object, the Ne estimates are labeled as `Estimated Ne`, and then confidence intervals are given using both a [parametric bootstrapping and jackknife method](https://www.datasciencecentral.com/resampling-methods-comparison/).
