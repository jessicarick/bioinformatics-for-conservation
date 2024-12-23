---
title: 13 - Detecting hybridization
type: docs
prev: docs/folder/12-landscape-genomics
next: docs/folder/14-molecular-forensics
---

## Detecting hybridization and introgression
Hybridization (mixing between two species) or admixture (mixing of genetics between two populations, not necessarily ones that are distinct species) can be a very complicated and complex topic, especially when viewed in a conservation and management context. Hybridization is not necessarily a bad thing in and of itself, and it can have both positive (e.g., increased genetic diversity) and negative (e.g., outbreeding depression or the extinction of a lineage) consequences for a population. In some situations, we may be interested in the presence or rate of hybridization itself, and in others we may want to check our dataset to make sure that we don't detect any individuals with mixed ancestry before proceeding with analyses that define individuals into distinct populations or groups.

Generally, there are two main types of methods used for detecting hybridization (1) clustering analyses, and (2) tree-based inferences. The former are often related to analyses used to detect population structure, such as the programs STRUCTURE or ADMIXTURE, which we discussed in [Week 9](https://jessicarick.github.io/bioinformatics-for-conservation/docs/folder/10-population-structure/). Admixed individuals often appear in these analyses as individuals that have partial genome assignment to multiple groups, and we can then use further analyses to sort these individuals into hybrid "classes" (e.g., F1, F2, F3, early-stage backcrosses, and late-stage backcrosses). With tree-based analyses, a hypothesized species tree phylogeny is usually used as the backbone for analyses, and statistics (e.g., [ABBA-BABA tests or D-statistics](https://avianhybrids.wordpress.com/2019/11/09/d-statistics-for-dummies-a-simple-test-for-introgression/)) are calculated to determine whether allele sharing between non-sister taxa is more than would be expected in the absence of admixture between them. Some of the programs and methods that employ these tree-based analyses include [HyDE](https://hybridization-detection.readthedocs.io/), [Dsuite](https://github.com/millanek/Dsuite), and [QuIBL](https://github.com/miriammiyagi/QuIBL). Which of these programs is best for your analyses will depend on what question you are most interested in, and what data you have available. For the purposes of our in-class exercises, we will be using the program [entropy](https://bitbucket.org/buerklelab/mixedploidy-entropy/src/master/), which is detailed below.

## Running entropy
We'll be using the program [entropy](https://bitbucket.org/buerklelab/mixedploidy-entropy/src/master/) to investigate the presence and identity of hybrids in the dataset from [Rosenthal et al. (2024)](https://doi.org/10.1002/ece3.11706), which sampled walleye (*Sander vitreus*) and sauger (*Sander candadensis*) from water bodies across Wyoming. For the purposes of our in-class exercises, I have already converted a VCF of these data (containing genotype likelihoods) to the input file format for entropy, which is a `.mpgl` file. This file is available in the `/xdisk/jrick/consbioinf/shared_data/week12_data/` directory and is named `entropy_in_sarwae_miss0.75_maf0.05_subsamp.mpgl`. This directory also contains a metadata file (which we'll be using for plotting the results in R), the script for plotting entropy results in R, and a file with starting q-value estimates `ldak_sarwae_2all.txt`, which is produced during the conversion of the VCF to mpgl format and makes the MCMC estimation process more efficient.

```sh
### FORMATTING INPUT
### on the UA HPC

module load R 

## convert VCF to mpgl input format
## requires installation of vcfR and MASS
## MASS is installed on the HPC, but vcfR may need to be
Rscript /xdisk/jrick/programs/entropy/auxfiles/inputdataformat.R sarwae_miss0.75_maf0.05_subsamp.vcf 

### RUNNING ENTROPY

## run entropy at k=2 (takes a while!! ~5-10min)
## if you copy this code, you'll have to remove the comments for it to run!
/xdisk/jrick/programs/entropy/entropy \
  -i entropy_in_sarwae_miss0.75_maf0.05_subsamp.mpgl \ # infile
	-n 2 \ # ploidy
	-l 3000 \ # mcmc chain length
	-b 1000 \ # burnin length
	-t 20 \ # only record output every 20 steps
	-k 2 \ # number of groups
	-Q 1 \ # estimate inter species ancestry (Q)
	-q ldak_sarwae_2all.txt \ # starting q value estimates (optional)
	-o entropy_out_sarwae.hdf5 # output file name

## post-process the results into text files
/xdisk/jrick/programs/entropy/estpost.entropy -p q -s 0 entropy_out_sarwae.hdf5 -o entropy_out_sarwae_k2_q.txt
/xdisk/jrick/programs/entropy/estpost.entropy -p Q -s 0 entropy_out_sarwae.hdf5 -o entropy_out_sarwae_k2_Q12.txt
```

Once we have these results, we can download them into our R-project `data/` folder and then use R to visualize the results!

```r
## wfsc 496b/596b, fall 2024
## week 12: hybridization
## post-processing and plotting entropy results

library(tidyverse)

inds <- read_table("data/sarwae_metadata.txt")

sarwae_q <- read_csv("data/entropy_out_sarwae_k2_q.txt") %>%
  mutate(pop = gsub("ind_[0-9]+_(pop_[0-9])","\\1",param),
         ind = gsub("q_(ind_[0-9]+)_pop_[0-9]+","\\1",param)) %>%
  pivot_wider(id_cols=ind, names_from=pop, values_from=median)

sarwae_Q <- read_csv("data/entropy_out_sarwae_k2_Q12.txt") %>%
  mutate(Qval = gsub("ind_[0-9]+_(anc_[0-9]-[0-9])","\\1",param),
         ind = gsub("Q_(ind_[0-9]+)_anc_[0-9]-[0-9]","\\1",param)) %>%
  pivot_wider(id_cols=ind, names_from=Qval, values_from=median)

# combine results together
sarwae_res <- inds %>%
  add_column(sarwae_q %>% dplyr::select(-ind),
             sarwae_Q %>% dplyr::select(-ind))

# first, plot an admixture barplot
# here, admixed individuals will have multiple colors
sarwae_res %>%
  pivot_longer(cols=starts_with("q_pop"),names_to="pop",values_to="qest") %>%
  ggplot(aes(x=ind, y=qest, fill=factor(pop))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_manual(values = c("#005f73","#f7b801"),
                    name="K",
                    labels=seq(1:2)) +
  facet_wrap(~sampling_loc, scales="free_x")

# now, a simple q-Q triangle plot
sarwae_res %>%
  ggplot(aes(x=q_pop_0, y=`Q_anc_2-1`)) +
  geom_point() +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1)

# now, more customization for the triangle plot (can play around with this!)
sarwae_res %>%
  ggplot(aes(x=q_pop_0, y=`Q_anc_2-1`)) +
  geom_segment(x=0,y=0,xend=0.5,yend=1,color="gray50",lty=3) +
  geom_segment(x=1,y=0,xend=0.5,yend=1,color="gray50",lty=3) +
  geom_segment(x=0,y=0,xend=1,yend=0,color="gray50",lty=3) +
  geom_point(aes(color=spp_id), size=4, alpha=0.5) +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  facet_wrap(~sampling_loc) +
  scale_color_manual(values = c("#005f73","#f7b801")) +
  xlab("Proportion walleye ancestry") +
  ylab("Interspecific ancestry (Q12)")

```

## Some notes
For the purposes of our in-class exercises, we only ran one MCMC chain and did not run it for very long. If you use entropy in the future, you will want to run multiple MCMC chains and assess their convergence to make sure that you're ending up with the same answer. You will also likely want to run the MCMC chain for much longer-- how long will depend on [an inspection of your MCMC traces](https://drvalle1.github.io/20_MCMC_convergence.html), but it will likely be for 50,000 steps or longer.
