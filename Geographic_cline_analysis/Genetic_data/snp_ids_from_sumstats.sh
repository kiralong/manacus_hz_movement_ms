#!/bin/bash

# Script to get input files for HZAR:
# 1. List of SNP IDS from Stacks sumstats file
# 2. List of scaffolds/chromosomes with more than 100 SNPs in them

# Populations sumstats
sms=populations.sumstats.tsv
# SNP Ids
ids=snp_ids_72k.tsv
# Scaffold with 100+ SNPs
scf=scaffolds_with_over100snps_list.tsv

# Get the SNP ids from the first four columns of the sumstats
cat $sms | grep -v '^#' | cut -f1-4 | sort -n -k1 -u > $ids

# From the SNP ids, tally the SNPs per-scaffold
cat $ids | cut -f2 | sort | uniq -c | awk '$1 >= 100 {print $2}' > $scf
