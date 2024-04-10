#!/bin/bash
#SBATCH -p aces
#SBATCH -N manacus_bwa_alignment
#SBATCH -t 168:00:00
#SBATCH -N 1
#SBATCH -n 8

# Load required modules for the cluster you are running on
module load gcc/7.2.0

# Declare paths to installed software
bwa_path=/path/to/bwa-0.7.17
samtools_path=/path/to/samtools-1.7
parallel_path=/path/to/parallel-20180322/bin

# Declare paths to input files and output locations
sample_names_path=/path/to/sample_info/barcodes_jan_2018.tsv
reads_path=/path/to/processed_samples
out_path=/path/to/output/dir/aligned_samples/ASM171598v3
database_path=/path/to/reference_genome/ASM171598v3/GCF_001715985.3_ASM171598v3_genomic.fna.gz

# Run through list of samples to align and sort all sample reads
cat $sample_names_path | cut -f 2 |
while read sample; do
	echo "$bwa_path/bwa mem -t 2 $database_path $reads_path/${sample}.1.fq.gz $reads_path/${sample}.2.fq.gz | \
		$samtools_path/samtools view -b -h | $samtools_path/samtools sort --threads 2 -o $out_path/${sample}.bam";
done | $parallel_path/parallel -j 4
