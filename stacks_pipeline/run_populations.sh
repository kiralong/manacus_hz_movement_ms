#!/bin/bash
#SBATCH -p catchenlab
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -J populations_job_name
#SBATCH -t 68:00:00

# Load required modules for computing cluster
module load gcc/7.2.0

# Declare paths to installed software, input, and output
stacks=/path/to/local/bin
popmap_path=/path/to/popmaps/full_manacus_min_6x_cov_no_dups_popmap.tsv
gstacks_out=/path/to/gstacks_runs/220308_gstacks_rm-pcr-dups_ALL
populations_output_path=/path/to/populations_runs
whitelist_path=/path/to/whitelists/whitelist.tsv

# Declare filtering parameters, as needed
mac=3
r=0.8
p=9
#maf=0.05

# Name and create output directory for the run
populations_output=$populations_output_path/populations_p${p}_r${r}_maf${maf}
mkdir -p $populations_output

# The populations command. Note the specific flags used change depending on what output you need, unused flags are commented out
cmd=(
    $stacks/populations
    --in-path $gstacks_out
    --out-path $populations_output
    --popmap $popmap_path
    --threads 6
    --min-samples-per-pop $r
    --min-mac $mac
#    --min-maf $maf
    --min-population $p
#    --hwe
#    --write-single-snp
#    --whitelist $whitelist_path
#    --fstats
#    --structure
#    --hzar
#    --genepop
#    --plink
    --vcf
    --ordered-export
#    --vcf-all
)

"${cmd[@]}"
