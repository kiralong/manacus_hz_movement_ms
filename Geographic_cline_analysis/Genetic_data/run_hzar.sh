#!/bin/bash
#SBATCH -p aces
#SBATCH -J hzar_Long_array
#SBATCH -t 168:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --array 1-101%6

# Bash script to run HZAR as an array job on a SLURM computing cluster

# Load version of R
module load R/4.0.1

# Declare paths and variables
thr=8
script_dir=/path/to/script/hzar/parallel_run
input_dir=$script_dir/Long
scaffold_list=$input_dir/scaffolds_with_over100snps_list.tsv
cd $input_dir

# Select a chromosome based on the Slurm array ID
chrom=$(cat $scaffold_list | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Make a new directory for that single run
out_dir=$input_dir/${chrom}
mkdir -p $out_dir
echo -e "Running HZAR for scaffold ${chrom}"

# HZAR script command
cmd=(
	$script_dir/hzarscript_loop_to_parallelize.R
	--workdir $out_dir
	--snpidfile $input_dir/snp_ids_72k.tsv
	--hzarfile $input_dir/populations.hzar.csv
	--chromosome $chrom
	--cores $thr
)

echo "${cmd[@]}"
"${cmd[@]}"
