# _Manacus_ hybrid zone movement manuscript

> Long _et al._ (2022) **Ongoing introgression of a secondary sexual plumage trait in a stable avian hybrid zone**. _bioRxiv_ <https://doi.org/10.1101/2023.03.30.535000>

Repository describing the bioinformatic code used for the analysis of population structure and movement of the _Manacus_ hybrid zone in western Panama.

## Stacks pipeline

Analysis of RADseq data using the [*Stacks*](https://catchenlab.life.illinois.edu/stacks/) software.

## Population structure analysis

PCA, admixture, intersp het, etc.

## Hybrid simulations

1. `make_diagnostic_snp_whitlist.py`: Generate a whitelist of a set of parental diagnostic SNPs filtered from a user-defined set of parameters. Takes a SUMSTATS and FSTATS output from *Stacks*.
2. `sim_hybrids_from_parents_vcf.py`: Take a VCF containing genotypes from two parental populations and simulate the genotypes across several hybrid crosses.

The original source for these RADseq-based hybrid simulations can be found in the following repository: <https://github.com/arcolon14/rad_hybrid_sims>.

## GGhybrid analysis

Stacks -> gghybrid format, hybrid index, etc.

## Geographic cline analysis (genetic data)

Stacks -> hzar format, hzar parallel script, hzar formatting and filtering

## Geographic cline analysis (phenotype data)

hzar phenotype scripts

## Data availability

Raw RAD-seq reads can be found on NCBI under BioProject [PRJNA893627](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA893627).

## Notes

### 2024-04-04

ARC created a local copy of the hybrid simulation scripts from <https://github.com/arcolon14/rad_hybrid_sims>.

### 2024-03-22

Started migration of the code from other repositories, e.g., <https://github.com/kiralong/RADseq_pipeline>, into a single, manuscript-specific repo.

## Authors

Kira M. Long
