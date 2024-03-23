#!/bin/bash
#SBATCH -p queque_name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J manacus_rad_process_ragtags
#SBATCH -t 168:00:00

module load gcc/7.2.0

raw=/path/to/raw/reads
out=/path/to/output/directory
bar=/path/to/barcodes.tsv

# Clean and demultiplex the samples.
cmd=(
    process_radtags
    -p $raw
    -o $out
    -b $bar
    --paired
    --clean
    --quality
    --rescue
    --renz-1 sbfI
    -i gzfastq
    --adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCAGAACAA   # P2 top from Hohenlohe2012
    --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT                # P1 bottom from Hohenlohe2012
    --adapter_mm 2
)

# Print command and run
echo "${cmd[@]}"
"${cmd[@]}"
