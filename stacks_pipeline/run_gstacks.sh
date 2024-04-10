#!/bin/bash
#SBATCH -p aces
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -J gstacks_manacus_ALL
#SBATCH -t 72:00:00

# Declare paths and variables
work_dir=/path/to/gstacks_runs
popmap=/path/to/popmaps/popmap_all_samples_pops_only_AN.tsv
aligned_samples=/path/to/aligned_samples/ASM171598v3/2021_dec_full_manacus_runs
out=$(date +${work_dir}/%y%m%d_gstacks_rm-pcr-dups_ALL)
threads=28

# Generate your time stamped gstacks output directory
mkdir -p $out

# Gstacks command to make catalog of genotypes and remove pcr duplicates
/path/to/programs/stacks-2.60/gstacks \
	-I $aligned_samples \
	-O $out \
	-M $popmap \
	-t $threads \
	--rm-pcr-duplicates \
	--details
