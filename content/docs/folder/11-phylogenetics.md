---
title: 11 - Phylogenetics
type: docs
prev: docs/folder/10-population-structure
next: docs/folder/12-landscape-genomics
---

## Building phylogenies
Phylogenies represent a hypothesis about the evolutionary history of a group of organisms, based on the data that we have on hand. These data can be morphological or genetic/genomic, and here we'll be discussing genomic-based phylogenetic reconstruction methods. There is a lot of important jargon related to phylogenies, which we discussed in class. To review this material, please watch the Zoom recording from Week 11; I am also including a picture of the white board information below:

![alt text](/content/img/whiteboard_phylogen_basics.jpeg | width=50%)
![alt text](/content/img/whiteboard_building_phylogenies.jpeg)

As with many of the other topics we've covered, there are many many options of ways to go about phylogenetic reconstruction. For better coverage of these, I'll refer you to a couple of nice review articles: X, X, and X. For most phylogenetic methods, we'll want to filter our variants more stringently than for other analyses, because the sites that will give us the most useful information for tree construction will be the sites that (1) have data for the majority of our taxa, and (2) have minor alleles that are shared by two or more individuals/taxa in our dataset. For the data we'll be working with here, I've filtered sites to only keep those with a minor allele frequency > 0.05, as well as only keeping sites with 10% or less missing data (i.e., `--max-missing 0.9`). 

For most tree building programs, we need to have our data in an alignment format-- generally, most programs take either a [`fasta`](https://zhanggroup.org/FASTA/), [`nexus`](https://plewis.github.io/nexus/), or [`phylip`](https://www.phylo.org/index.php/help/phylip) format. I've already converted the filtered VCF we'll be using into a `.phy` file; if you need to do this in the future on your own, the easiest way is to use the `vcf2phylip.py` script from [this GitHub repository](https://github.com/edgardomortiz/vcf2phylip) (if you need this for your Project 2, I have also included it in the `project2_data/` directory in our `shared_data/` folder). This script can be used as follows:

```sh
python vcf2phylip.py -i gila_variants_outgroup_filtered.vcf 
```

This will output a phylip file with the same name as the input VCF file. There are other options for customizing the conversion process (including filtering again for missing data or minor allele frequency), which you can explore using `python vcf2phylip.py --help`.

For our purposes, we'll be working with [IQtree](http://www.iqtree.org/) due to its simplicity and flexibility, which means that it's really easy to do very basic tree construction, but there are also a ton of options for doing much fancier things as well. We'll just be building a simple phylogeny from the `gila_variants_filtered_outgroup.phy` data from the `week10_data/` directory, and then using `R` to visualize our results. These data contain 71 sequences, which includes 70 ingroup taxa and 1 outgroup individual ("SRR15431420"). To run IQtree using the basic settings, we'll use the following code:

```sh
module load iqtree2

iqtree2 -s gila_variants_filtered_outgroup.phy -st DNA --prefix gila_tree
```
There will be a lot of output to the screen about what is occurring as IQtree runs, which will also be saved in some of the output files. If we wanted to additionally run 1000 bootstraps to quantify our confidence in the tree, we can add a `-B 1000` flag to the `iqtree2` code. For our purposes in-class, we won't be doing this because the bootstrap analyses take a while to run.

From here, we'll want to download the `.treefile` output file from IQtree (this is our best estimate for the phylogeny for these taxa), as well as the metadata file (`gila_selected_inds.txt` from the `week10_data/` directory). Just like previously, we'll want to download these files to the `data/` directory in our R project, using either `scp` or the interactive OOD online interface. Once we have these files downloaded, we can use them to visualize our phylogeny in R. To do so, we'll be using three new-to-us packages, `ape`, `phytools`, and `ggtree`; these may take a bit of time to install if you do not already have them.

```r
#install.packages("ape")
#install.packages("phytools")
#install.packages("ggtree")
library(ape)
library(phytools)
library(ggtree)
library(tidyverse)

# load sequence metadata
ind_dat <- read_csv("data/gila_selected_inds.txt") %>%
  select(Organism, Ecotype, Run, `Sample Name`) %>%
  mutate(Pop = gsub("[0-9]+([A-Z]+)[0-9]+","\\1",`Sample Name`))

# load tree
tre <- read.tree("data/gila_tree.treefile")

# root tree using the outgroup
tre <- root(tre,"SRR15431420",
            resolve.root=TRUE)

# plot simple tree
plot.phylo(tre, show.tip.label=FALSE)

# plot using ggtree
ggtree(tre)

# match metadata to tree
tre_metadat <- tibble(bam=tre$tip.label) %>%
  mutate(ind=bam) %>%
  left_join(ind_dat, by=c("ind"="Run")) 

p <- ggtree(tre, layout="equal_angle") %<+% tre_metadat
p + 
  geom_tippoint(aes(color=Organism)) +
  geom_tiplab(aes(label=Pop), color="black")
```
