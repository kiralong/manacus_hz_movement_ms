#!/usr/bin/env python3
import sys, os, argparse, gzip, random, datetime, warnings

# (c) 2023 Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]
DESC = '''Identify a set of diagnostic SNPs between two parental
populations based on the FSTATS and SUMSTATS outputs of STACKS.'''

def parse_args():
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-s', '--sumstats',
                   required=True,
                   help='(str) Path to the SUMSTATS file from POPULATIONS.')
    p.add_argument('-f', '--fstats', required=True,
                   help='(str) Path to the FST_Y-Z file from POPULATIONS.')
    p.add_argument('-o', '--outdir',
                   required=False, default='.',
                   help='Output directory.  [str]')
    p.add_argument('-n', '--n-sites', required=False, default=1000, type=int,
                   help='(int) Number of individuals simulated per population. [default=1000]')
    p.add_argument('-e', '--hwe',
                   required=False, action='store_true', default=False,
                   help='(bool) Discard sites out of HWE, applied per population [default=False]')
    p.add_argument('-a', '--hwe-alpha',
                   required=False, default=0.05, type=float,
                   help='(float) p-value cutoff to clasify a site as out of HWE [default=0.05]')
    p.add_argument('-p', '--private',
                   required=False, action='store_true', default=False,
                   help='(bool) Keep only sites containing private alleles in the parental populations [default=False]')
    p.add_argument('-d', '--min-fst',
                   required=False, default=0.25, type=float,
                   help='(float) Min Fst value used to classify a site as divergent between parental populations [default=0.25]')
    p.add_argument('-m', '--min-maf', required=False, default=0.05, type=float,
                   help='(float) Minimum allele frequency to retain a parental allele [default=0.05]')
    p.add_argument('-r', '--write-random-snp',
                   required=False, action='store_true', default=False,
                   help='(bool) Export only one random SNP per locus [default=False]')
    # Check input arguments
    args = p.parse_args()
    args.outd = args.outdir.rstrip('/')
    if not os.path.exists(args.sumstats):
        sys.exit(f"Error: '{args.sumstats}' not found.")
    if not os.path.exists(args.fstats):
        sys.exit(f"Error: '{args.fstats}' not found.")
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if args.n_sites <= 0:
        sys.exit(f"Error: 'number-sites' ({args.n_sites}) must be non-zero positive integer.")
    return args


def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'


def initialize_log(args, log=sys.stdout):
    '''Initialize the log file'''
    print(f'''{PROG} started on {now()}\n
Input parameters:
    SUMSTATS file:  {args.sumstats}
    FSTATS file:    {args.fstats}
    Output dir:     {args.outdir}
    N sites:        {args.n_sites:,}
    Min MAF:        {args.min_maf:0.06g}
    Min Fst:        {args.min_fst:0.06g}''',
    file=log, flush=True)
    if args.private:
        print('    Retaining private alleles', file=log, flush=True)
    if args.hwe:
        print(f'    Discarding sites out of HWE (p<{args.hwe_alpha:0.06g})',
              file=log, flush=True)
    if args.write_random_snp:
        print('    Exporting one random SNP per RAD locus', 
              file=log, flush=True)


class FiltParams:
    '''Filter parameters'''
    def __init__(self, n_pops=2, hwe=False, hwe_alpha=0.05, min_maf=0.05, min_fst=0.25, private=False):
        # Check inputs
        assert type(n_pops) in [int, float]
        assert n_pops >= 1, f'n_pops ({n_pops}) must be greater or equal than 1.'
        # HWE
        assert type(hwe) is bool
        assert type(hwe_alpha) is float, f'Error: hwe_alpha ({hwe_alpha}) not a floating point number.'
        assert 0.0 < hwe_alpha <= 1.0, f'HWE alpha ({hwe_alpha}) must be between 0 and 1.'
        # MAF
        assert type(min_maf) is float, f'Error: min_maf ({min_maf}) not a floating point number.'
        assert 0.0 <= min_maf <= 0.5, f'Min MAF ({min_maf}) must be between 0.0 and 0.5.'
        # Fst
        assert type(min_fst) is float, f'Error: min_fst ({min_fst}) not a floating point number.'
        assert 0.0 <= min_fst <= 1.0, f'min_fst ({min_fst}) must be between 0 and 1.'
        # Private alleles
        assert type(private) is bool
        # Assign parameters
        self.pops  = n_pops
        self.hwe   = hwe
        self.alpha = hwe_alpha
        self.maf   = min_maf
        self.fst   = min_fst
        self.priv  = private
    def filter_minor_allele(self, p):
        remove = False
        # Filter by MAF
        if p >= 0.5:
            if 0 < 1-p < self.maf:
                remove = True
        else:
            if 0 < p < self.maf:
                remove = True
        return remove


def parse_filter_params(prog_opt):
    '''Generate the class with the filtering options'''
    filt_opts = FiltParams(
        n_pops=2,  # Harcoded, we ALWAYS want two parental pops
        hwe=prog_opt.hwe,
        hwe_alpha=prog_opt.hwe_alpha,
        min_maf=prog_opt.min_maf,
        min_fst=prog_opt.min_fst,
        private=prog_opt.private)
    return filt_opts


class SumstatsSite:
    '''Site extracted from the SUMSTATS file'''
    def __init__(self, locus_id, snp_col, chrom, bp, popid, p, hwe, private):
        self.locid  = locus_id
        self.snpcol = snp_col
        self.chrom  = chrom
        self.bp     = bp
        self.popid  = popid
        self.p      = p
        self.hwe    = hwe
        self.priv   = private
    def __str__(self):
        return f'{self.locid}\t{self.snpcol}\t{self.chrom}\t{self.bp}\t{self.popid}\t{self.p}\t{self.hwe}\t{self.priv}'


def print_sumstats_tally(sites_dict, log):
    '''Print a tally of filtered loci from SUMSTATS'''
    # Tally total kept sites
    nlocs = len(sites_dict)
    nsnps = 0
    for loc in sites_dict:
        for site in sites_dict[loc]:
            nsnps += 1
    print('\nAfter filtering the SUMSTASTS, the program kept:', file=log)
    print(f'    Total loci: {nlocs:,}', file=log)
    print(f'    Total SNPs: {nsnps:,}', file=log, flush=True)


def parse_sumstats(sumstats_f, filt_opts, log):
    '''Parse FSTATS file'''
    assert isinstance(filt_opts, FiltParams)
    # Print to log
    sum_f = sumstats_f.split('/')[-1]
    print(f'\nProcessing {sum_f}...', file=log, flush=True)
    # Output
    sites_dict = dict()
    prev_snp = None
    pops_in_site = list()
    # Tally
    record = 0
    loci = set()
    snps = set()
    # Parse the sumstats
    with open(sumstats_f, 'r') as fh:
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            # Parse the rows
            fields  = line.strip('\n').split('\t')
            locid   = int(fields[0])      # Locus ID
            snpcol  = int(fields[3])      # SNP column
            chrom   = fields[1]           # Chromosome ID
            bp      = int(fields[2])      # Basepair
            popid   = fields[4]           # Population ID
            p_freq  = float(fields[8])    # P allele freq
            hwe_p   = float(fields[19])   # HWE p-value
            priv    = int(fields[20])     # Private allele status
            snp     = f'{locid}_{snpcol}'  # SNP ID
            record += 1
            loci.add(locid)
            snps.add(snp)
            # Filter by maf/mac
            rm_by_allele = filt_opts.filter_minor_allele(p_freq)
            if rm_by_allele:
                continue
            # Filter by HWE
            if filt_opts.hwe:
                if hwe_p == -1.0:
                    # The HWE was not used in populations
                    warnings.warn('HWE P-value is \'-1.00000\'. The --hwe was not added to populations. Skipping HWE filtering.')
                    filt_opts.hwe = False
                    hwe_p = 1.0
                if hwe_p < filt_opts.alpha:
                    continue
            # Filter by private alleles
            if filt_opts.priv:
                if priv != 1:
                    continue
            # Process remaining sites
            site = SumstatsSite(locid, snpcol, chrom, bp, popid, p_freq, hwe_p, priv)
            # For the first snp
            if prev_snp is None:
                prev_snp = snp
            # If looking at the same SNP in a different pop just add to the list
            if snp == prev_snp:
                pops_in_site.append(site)
            # When looking at a different SNP
            else:
                prev_loc = int(prev_snp.split('_')[0])
                prev_col = int(prev_snp.split('_')[1])
                if len(pops_in_site) >= filt_opts.pops:
                    sites_dict.setdefault(prev_loc, {})
                    sites_dict[prev_loc][prev_col] = pops_in_site
                pops_in_site = list()
                prev_snp = None
                pops_in_site.append(site)
            prev_snp = snp
        prev_loc = int(prev_snp.split('_')[0])
        prev_col = int(prev_snp.split('_')[1])
        if len(pops_in_site) >= filt_opts.pops:
            sites_dict.setdefault(prev_loc, {})
            sites_dict[prev_loc][prev_col] = pops_in_site
    # Print to log
    print(f'    Read {record:,} total records from input SUMSTATS.', file=log)
    print(f'    Processed {len(loci):,} total input loci, composed of {len(snps):,} total SNPs.', file=log)
    print_sumstats_tally(sites_dict, log)
    return sites_dict


class FstatsSite:
    '''Site extracted from FSTATS file'''
    def __init__(self, locus_id, snp_col, chrom, bp, pop_a, pop_b, fst, fst_pval):
        self.locid  = locus_id
        self.snpcol = snp_col
        self.chrom  = chrom
        self.bp     = bp
        self.popa   = pop_a
        self.popb   = pop_b
        self.fst    = fst
        self.fstp   = fst_pval
    def __str__(self):
        return f'{self.locid}\t{self.snpcol}\t{self.chrom}\t{self.bp}\t{self.popa}\t{self.popb}\t{self.fst}'


def print_fstats_tally(sites_dict, log):
    '''Print a tally of filtered loci from FSTATS'''
    # Tally total kept sites
    nlocs = len(sites_dict)
    nsnps = 0
    for loc in sites_dict:
        for site in sites_dict[loc]:
            nsnps += 1
    print('\nAfter filtering the FSTATS, the program kept:', file=log)
    print(f'    Total loci: {nlocs:,}', file=log)
    print(f'    Total SNPs: {nsnps:,}', file=log, flush=True)


def parse_fstats(fstats_f, filt_opts, log):
    '''Parse FSTATS file'''
    assert isinstance(filt_opts, FiltParams)
    # Print to log
    fst_f = fstats_f.split('/')[-1]
    print(f'\nProcessing {fst_f}...', file=log, flush=True)
    # Output
    sites_dict = dict()
    # Tally
    records = 0
    loci = set()
    snps = set()
    # Parse the file
    with open(fstats_f, 'r') as fh:
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            # Parse the rows
            records += 1
            fields = line.strip('\n').split('\t')
            locid  = int(fields[0])       # Locus ID
            snpcol = int(fields[5])       # SNP column
            chrom  = fields[3]            # Chromosome
            bp     = int(fields[4])       # Basepair
            pop_a  = fields[1]            # Pop A ID
            pop_b  = fields[2]            # Pop B ID
            fst    = float(fields[11])    # AMOVA Fst
            fpval  = float(fields[13])    # AMOVA Fst Pval
            snp    = f'{locid}_{snpcol}'  # SNP ID
            loci.add(locid)
            snps.add(snp)
            # Filter by Fst
            if fst < filt_opts.fst:
                continue
            # TODO: Apply a Fst-pval filter???
            site = FstatsSite(locid, snpcol, chrom, bp, pop_a, pop_b, fst, fpval)
            # Add to the site dictionary
            sites_dict.setdefault(locid, dict())
            sites_dict[locid].setdefault(snpcol, dict())
            sites_dict[locid][snpcol] = site
    # Print to log
    print(f'    Read {records:,} total records from input FSTATS.', file=log)
    print(f'    Processed {len(loci):,} total input loci, composed of {len(snps):,} total SNPs.', file=log)
    print_fstats_tally(sites_dict, log)
    return sites_dict


def sample_kept_sites(kept_snps, n_sites=1000, outd='.', single_snp=False, log=sys.stdout):
    '''Sample the kept sites and save to the whistlist file'''
    n_loci = list(kept_snps.keys())
    export = 0
    if len(n_loci) == 0:
        return
    with open(f'{outd}/snp_whitelist.tsv', 'w') as fh:
        # Exporting single snp per-locus
        if single_snp:
            print('\nExporting sites (random single SNP per locus)...', file=log, flush=True)
            if n_sites > len(n_loci):
                print(f'    Warning: more sites chosen for export ({n_sites:,}) than loci available post filtering ({len(n_loci):,}). Exporting {len(n_loci):,} total sites.', file=log, flush=True)
                n_sites = len(n_loci)
            sampled_loci = random.sample(n_loci, n_sites)
            for locus in sorted(sampled_loci):
                snp_cols = kept_snps[locus]
                assert len(snp_cols) > 0, f'{locus} {snp_cols}'
                col = random.choice(snp_cols)
                fh.write(f'{locus}\t{col}\n')
            print(f'    Exported {n_sites:,} total sites to the whitelist.', file=log, flush=True)
        # When exporting multiple SNPs per locus:
        else:
            print('\nExporting sites (more than one variant sites can be exported per locus)...', file=log, flush=True)
            # Make a new list of SNP ids, regardless of locus
            snps = list()
            for loc in sorted(kept_snps):
                for col in sorted(kept_snps[loc]):
                    snp_id = (loc, col)
                    snps.append(snp_id)
            # Adjust if the user asks for more sites than available
            if n_sites > len(snps):
                print(f'    Warning: more sites chosen for export ({n_sites:,}) than sites available post filtering ({len(snps):,}). Exporting {len(snps):,} total sites.', file=log, flush=True)
                n_sites = len(snps)
            sampled_snps = random.sample(snps, n_sites)
            for snp in sorted(sampled_snps):
                fh.write(f'{snp[0]}\t{snp[1]}\n')
            print(f'    Exported {len(sampled_snps):,} total sites to the whitelist.', file=log, flush=True)
    fh.close()


def compare_sites(sumstats_sites, fst_sites, log):
    '''Find the overlap in the SNPs that are retained across both datasets'''
    print('\nMerging retained sites for SUMSTATS and FSTATS...', file=log, flush=True)
    sums_loc_set = set(sumstats_sites.keys())
    fsts_loc_set = set(fst_sites.keys())
    kept_loci = sorted(list(sums_loc_set.intersection(fsts_loc_set)))
    nsnps = 0
    kept_snps = dict()
    for locus in kept_loci:
        fst_loc_snps = set(fst_sites[locus].keys())
        sum_loc_snps = set(sumstats_sites[locus])
        olap_snps = sorted(list(sum_loc_snps.intersection(fst_loc_snps)))
        if len(olap_snps) == 0:
            continue
        nsnps += len(olap_snps)
        kept_snps[locus] = olap_snps
    print(f'    Found {len(kept_snps):,} loci (composed of {nsnps:,} SNPs) shared across both dataset.', file=log, flush=True)
    return kept_snps


def main():
    args = parse_args()
    # Initialize log file
    log_f = open(f'{args.outdir}/diagnostic_parental_wl.log', 'w')
    initialize_log(args, log_f)
    # Prepare the filtering options
    filt_opts = parse_filter_params(args)
    # Process SUMSTATS
    sumstats_sites = parse_sumstats(args.sumstats, filt_opts, log_f)
    # Process FSTATS
    fstats_sites = parse_fstats(args.fstats, filt_opts, log_f)
    # Compare the two loci
    kept_snps = compare_sites(sumstats_sites, fstats_sites, log_f)
    # Generate the whitelist file
    sample_kept_sites(kept_snps, args.n_sites, args.outdir,
                      args.write_random_snp, log_f)

    # Close outputs
    log_f.write(f'\n{PROG} finished on {now()}\n')
    log_f.close()

# Run Code
if __name__ == '__main__':
    main()
