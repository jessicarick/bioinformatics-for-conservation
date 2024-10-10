## Estimating relatedness, kinship, and parentage
All three of these tasks are somewhat related, but can be approached in several different ways. Essentially, they all entail quantifying the relationships between individuals in a dataset, and then in the case of kinship and parentage, then binning those relationships into parent-offpring, full sibling, half sibling, etc. according to theory on how related individuals should be for each of these categories. The most recent methods for these sorts of estimation methods use Bayesian or maximum likelihood inference to quantify our uncertainty in a given estimate of kinship or a given relationship between individuals. 

For our class, we're going to approach this two different ways. First, we're going to estimate relatedness among all individuals in a dataset using [NgsRelate](https://github.com/ANGSD/NgsRelate), which is part of the ANGSD ecosystem of programs. NgsRelate will actually co-estimate both relateness and inbreeding coefficients, so we can efficiently get estimates for both of these measures using the same program. We can start analyses in NgsRelate from either `.bam` or `.vcf` files (or from a `.glf` file output from ANGSD), depending on our analysis pipeline.

For the simplest estimation of relatedness among individuals, we can use the following code:

```sh

```
