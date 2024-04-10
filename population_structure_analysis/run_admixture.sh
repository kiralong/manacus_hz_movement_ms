#!/bin/bash
#SBATCH -p aces
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -J Long_10K_admixture
#SBATCH -t 168:00:00

# Set variables and paths
threads=4
working_dir=/path/to/work/dir/admixture/Long_10k
plink_basename=$working_dir/populations.corrected.plink

# Runs 9 k clusters in admixture and generates output directories for each k
for k in {1..9}
do
    out_dir=$working_dir/K${k}
    mkdir -p $out_dir
    cd $out_dir
    admixture --cv -j$threads ${plink_basename}.bed $k > $out_dir/${k}.log
    paste ${plink_basename}.fam $out_dir/populations.corrected.plink.${k}.Q | tr ' ' '\t' | cut -f 1,2,7- > $out>
done
