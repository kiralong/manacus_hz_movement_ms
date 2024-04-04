#!/usr/bin/env python3
# (c) 2023 Angel G. Rivera-Colon
import sys, os, gzip, argparse, datetime

PROG = sys.argv[0].split('/')[-1]

#
# Command line options
#
def parse_args():
    desc = '''Format a VCF into a genotype matrix \
    compatible with the INTROGRESS R package.'''
    p = argparse.ArgumentParser(prog=PROG, description=desc)
    p.add_argument('-v', '--vcf',
                   required=True, help='(str) Path to the input VCF.')
    p.add_argument('-p', '--popmap',
                   required=True, help='(str) Path to input POPMAP')
    p.add_argument('-o', '--outd',
                   required=False, default='.',
                   help='(str) Path to output directory [default=./]')
    # Check input arguments
    args = p.parse_args()
    args.outd = args.outd.rstrip('/')
    if not os.path.exists(args.vcf):
        sys.exit(f"Error: '{args.vcf}' not found.")
    if not os.path.exists(args.popmap):
        sys.exit(f"Error: '{args.popmap}' not found.")
    if not os.path.exists(args.outd):
        sys.exit(f"Error: '{args.outd}' not found.")
    return args

def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def parse_popmap(popmap_f):
    '''Parse the population map and get population asignment of the individuals.'''
    print('Parsing POPMAP...')
    # List of populations
    pops = list()
    # Key value pair of individual population asignment
    population_map = dict()
    with open(popmap_f) as fh:
        for line in fh:
            fields = line.strip('\n').split('\t')
            indv = fields[0]
            pop = fields[1]
            # Add to dictionary
            population_map[indv] = pop
            # Tally populations
            if pop not in pops:
                pops.append(pop)
    print(f'    Read {len(population_map):,} individuals across {len(pops):,} populations.\n')
    return population_map

def parse_vcf(vcf_f, popmap, outdir):
    '''Parse the VCF to extract the genotypes and export output files.'''
    print('Parsing VCF...')
    # Prepare output files
    loci_f = open(f'{outdir}/loci_tbl.csv', 'w')
    loci_f.write('Locus,type\n')
    geno_f = open(f'{outdir}/genotypes_tbl.csv', 'w')
    vcf_rec = 0
    # Open VCF
    with gzip.open(vcf_f, 'rt') if vcf_f.endswith('.gz') else open(vcf_f) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                vcf_samples = line.strip('\n').split('\t')[9:]
                if len(vcf_samples) != len(popmap):
                    sys.exit(f"Error: Number of samples in popmap ({len(popmap)}) and input VCF ({len(vcf_samples)}) do not match.")
                print(f'    Processing genotypes from {len(vcf_samples):,} input samples.')
                # Preparet the pop and sample header for genotype tbl
                pop_row = list()
                for sam in vcf_samples:
                    pop = popmap.get(sam, None)
                    if pop is None:
                        sys.exit(f"Error: POPMAP sample {sam} not in VCF samples.")
                    pop_row.append(pop)
                pop_row = ','.join(pop_row) + '\n'
                sam_row = ','.join(vcf_samples) + '\n'
                geno_f.write(pop_row)
                geno_f.write(sam_row)
            else:
                # Process the VCF records
                vcf_rec += 1
                fields = line.strip('\n').split('\t')
                # Extract the locus info for loci.tbl
                info = fields[2]
                loc = info.split(':')[0]
                col = info.split(':')[1]
                loci_f.write(f'{loc}_{col},C\n')
                # Process the genotypes
                genotypes = fields[9:]
                geno_row = list()
                for geno in genotypes:
                    # Format missing genotypes
                    if geno in ['./.', '.|.', '.']:
                        geno = 'NA/NA'
                    geno_row.append(geno)
                geno_row = ','.join(geno_row) + '\n'
                geno_f.write(geno_row)
    print(f'    Formated genotypes for {vcf_rec:,} records.')
    geno_f.close()
    loci_f.close()

def parse_gghybrid_csv(gghybrid_f, outdir='.'):
    '''Parse the GGHYBRID CVS to extract genotypes.'''
    # Outputs
    ind_pops = dict()   # Individual-Population key-value pairs
    ind_genos = dict()  # Genotypes of a single individual
    loci = list()       # Locus IDs
    # Open input
    fh = open(gghybrid_f, 'r')
    if gghybrid_f.endswith('gz'):
        fh = gzip.open(gghybrid_f, 'rt')
    for i, line in enumerate(fh):
        if line.startswith('#'):
            continue
        if line.startswith('INDLABEL'):
            # First data line is: "INDLABEL,POPID,LOCUS_1,LOCUS_2,...,LOCUS_N"
            loci = line.strip('\n').split(',')[2:]
            print(f'Read {len(loci):,} input loci.')
        else:
            fields = line.strip('\n').split(',')
            indv = fields[0]
            pop  = fields[1]
            # Add to the dictionary of individuals+pops
            ind_pops[indv] = pop
            # Get the genotypes for the individual
            if indv not in ind_genos:
                # The first time you see the sample you have to initialize the
                # data structure and add the first allele for each locus
                ind_genos.setdefault(indv, [[], []])
                genos = fields[2:]
                assert len(genos) == len(loci)
                ind_genos[indv][0] = genos
            else:
                # The sample is already present, so just add the other allele
                genos = fields[2:]
                assert len(genos) == len(loci)
                ind_genos[indv][1] = genos
    print(f'Read data for {len(ind_pops.keys()):,} individuals across {len(set(ind_pops.values())):,} populations.')
    fh.close()
    # Generate Outputs
    # 1. Loci table
    loci_f = f'{outdir}/loci_tbl.csv'
    with open(loci_f, 'w') as fh:
        fh.write('Locus,type\n')
        for loc in loci:
            row = f'{loc},C\n' # C because the markers are "co-dominatn" for INTROGRESS
            fh.write(row)
        print(f'Saved loci table to:\n\t{loci_f}')
        fh.close()
    # Genotype matrix
    geno_f = f'{outdir}/genotypes_tbl.csv'
    with open(geno_f, 'w') as fh:
        # Add the headers
        pop_row = list()
        ind_row = list()
        for pair in ind_pops.items():
            pop_row.append(pair[1])
            ind_row.append(pair[0])
        pop_row = ','.join(pop_row)
        ind_row = ','.join(ind_row)
        fh.write(f'{pop_row}\n')
        fh.write(f'{ind_row}\n')
        # Add the genotypes
        for l, loc in enumerate(loci):
            locus_genos = list()
            for ind in ind_pops:
                ind_geno = ind_genos[ind]
                a1 = ind_geno[0][l]
                a2 = ind_geno[1][l]
                locus_genos.append(f'{a1}/{a2}')
            genos_row = ','.join(locus_genos)
            fh.write(f'{genos_row}\n')
        print(f'Saved genotypes table to:\n\t{geno_f}')
        fh.close()



def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    population_map = parse_popmap(args.popmap)
    genotypes = parse_vcf(args.vcf, population_map, args.outd)
    print(f'\n{PROG} finished on {now()}')

# Run Code
if __name__ == '__main__':
    main()
