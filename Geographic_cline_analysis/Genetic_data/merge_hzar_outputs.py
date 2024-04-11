#!/usr/bin/env python3
import os, sys, argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-c','--chromosome-list',required=True,help='Chromosome list')
    p.add_argument('-a','--hzar-run-a',required=True,help='Path to HZAR run A')
    p.add_argument('-b','--hzar-run-b',required=True,help='Path to HZAR run B')
    args = p.parse_args()
    if not os.path.exists(args.chromosome_list):
        sys.exit("Error: Chromosome list not found.")
    if not os.path.exists(args.hzar_run_a):
        sys.exit("Error: Path to HZAR run A not found.")
    if not os.path.exists(args.hzar_run_b):
        sys.exit("Error: Path to HZAR run B not found.")
    return args

def load_chromosomes(chrom_f):
    '''Make a list of all the chromosomes to process all the runs individuallly'''
    chromosomes = list()
    for lines in open(chrom_f, 'r'):
        chrom = lines.strip('\n')
        chromosomes.append(chrom)
    return chromosomes

def parse_hzar_output_file(hzar_f):
    '''Process a single hzar output table'''
    hzar_dict = dict()
    for i, lines in enumerate(open(hzar_f, 'r')):
        if lines.startswith('chromosome'):
            continue
        fields = lines.strip('\n').split('\t')
        snp = fields[1]
        model = fields[3]
        if model == 'nullModel':
            continue
        hzar_dict[snp] = fields
        # if i>25:
        #     break
    return hzar_dict

def get_unique_snps(a_snps, b_snps):
    '''Get the unique snps across both runs'''
    snps_list = list()
    for snp in a_snps:
        snps_list.append(snp)
    for snp in b_snps:
        snps_list.append(snp)
    snp_set = set(snps_list)
    return snp_set

def merge_two_cline_dicts(dict_a, dict_b):
    '''Merge the two HZAR runs by the matching SNPs'''
    # Get the individual SNPs
    snps = get_unique_snps(dict_a.keys(), dict_b.keys())
    # Get the pair of clines for each snps
    cline_pairs = dict()
    for snp in sorted(snps):
        # Get cline A for the snp
        cline_a = dict_a.get(snp, None)
        if cline_a is None:
            continue
        # Get cline B for the snp
        cline_b = dict_b.get(snp, None)
        if cline_b is None:
            continue
        pair = (cline_a, cline_b)
        cline_pairs[snp] = pair
    return cline_pairs

def prepare_output_table(run_dir):
    '''Initialize the new output table'''
    tsv_f = f'{run_dir}/snp_cline_parameters.filtered.tsv'
    fh = open(tsv_f, 'w')
    header = 'chromosome\tsnp_id\tbase_pair\tbest_model\tcenter\twidth\tpMin\tpMax\tcenter_low\tcenter_med\tcenter_high\twidth_low\twidth_med\twidth_high\n'
    fh.write(header)
    return fh

def save_cline_to_file(cline_list, fh):
    assert len(cline_list) == 14
    cline_str = '\t'.join(cline_list)
    fh.write(f'{cline_str}\n')

def process_single_chromosome(chromosome, run_a_dir, run_b_dir, a_fh, b_fh):
    '''Process the output for a single chromosome across the two parallel HZAR runs'''
    # Process the first run, aka A
    table_a = f'{run_a_dir}/{chromosome}/snp_cline_parameters.tsv'
    a_dict = dict()
    if os.path.exists(table_a):
        a_dict = parse_hzar_output_file(table_a)
    # Process the second run, aka B
    table_b = f'{run_b_dir}/{chromosome}/snp_cline_parameters.tsv'
    b_dict = dict()
    if os.path.exists(table_b):
        b_dict = parse_hzar_output_file(table_b)
    # Merge the dicts
    cline_pairs = merge_two_cline_dicts(a_dict, b_dict)
    # Save into the corresponding files
    for snp in sorted(cline_pairs):
        cline_pair = cline_pairs[snp]
        # Save the cline A
        save_cline_to_file(cline_pair[0], a_fh)
        # Save the cline B
        save_cline_to_file(cline_pair[1], b_fh)

def process_all_data(chroms_f, run_a_dir, run_b_dir):
    '''Process all the data'''
    # Load chromosomes
    chromosomes = load_chromosomes(chroms_f)
    # Create the two output files
    a_fh = prepare_output_table(run_a_dir)
    b_fh = prepare_output_table(run_b_dir)
    # Loop over the chromosomes and process
    for chromosome in sorted(chromosomes):
        process_single_chromosome(chromosome, run_a_dir, run_b_dir, a_fh, b_fh)
    a_fh.close()
    b_fh.close()

def main():
    args = parse_args()
    process_all_data(args.chromosome_list, args.hzar_run_a, args.hzar_run_b)

if __name__ == "__main__":
    main()
