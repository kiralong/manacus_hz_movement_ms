#!/usr/bin/env python3
import sys, os, argparse, random, warnings, datetime

# (c) 2023 Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]

#
# Command line options
#
def parse_args():
    desc = '''Filter a SUMSTATS file from the STACKS POPULATIONS program (populations.sumstats.tsv) \
        to generate a whitelist containing the STACKS catalog IDs of a subset of selected SNPs. The \
        sites in the SUMSTATS table can be filtered to remove sites out of Hardy-Heinberg equilibrium, \
        under a given allele frequency/count cutoff, and under a given number of specified populations/\
        samples. The user can specify the number of sites to be retained after filtering. The final \
        whitelist can export one SNP per locus, equivalent to the \'write-random-snp\' option in POPULATIONS.'''
    p = argparse.ArgumentParser(prog=PROG, description=desc)
    p.add_argument('-s', '--sumstats',
                   required=True, help='(str) Path to the populations sumstats TSV file')
    p.add_argument('-o', '--outd',
                   required=False, default='.',
                   help='(str) Path to output directory [default=./]')
    p.add_argument('-p', '--min-pops',
                   required=False, default=1, type=int,
                   help='(int) Minimum number of populations required to retain a site [default=1]')
    p.add_argument('-r', '--min-perc-indv',
                   required=False, default=None,
                   help='(float) Minimum percentage of individuals in a population required to retain a site [default=None]')
    p.add_argument('-i', '--min-num-indv',
                   required=False, default=None,
                   help='(int) Minimum number of individuals in a population required to retain a site. Cannot be used alongside \'min-perc-indv\' [default=None]')
    p.add_argument('-n', '--number-sites',
                   required=False, default=1000, type=int,
                   help='(int) Number of sites to export [default=1000]')
    p.add_argument('-e', '--hwe',
                   required=False, action='store_true', default=False,
                   help='(bool) Retain only sites in in HWE, applied per population [default=False]')
    p.add_argument('-a', '--hwe-alpha',
                   required=False, default=0.05, type=float,
                   help='(float) p-value cutoff to clasify a site as out of HWE [default=0.05]')
    p.add_argument('-f', '--min-maf',
                   required=False, default=None,
                   help='(float) Minimum allele frequency cutoff to retain a site, applied per-population. Cannot be used alongside \'min-mac\'. [default=None]')
    p.add_argument('-c', '--min-mac',
                   required=False, default=None,
                   help='(int) Minimum allele count cutoff to retain a site, applied per-population. Cannot be used alongside \'min-mac\'. [default=None]')
    p.add_argument('-t', '--max-obs-het',
                   required=False, default=1.0, type=float,
                   help='(float) Max observed heterozygosity cutoff to retain a site, applied per-population [default=1.0, no cutoff]')
    p.add_argument('-g', '--write-random-snp',
                   required=False, action='store_true', default=False,
                   help='(bool) Export only one random SNP per locus [default=False]')
        
    # Check input arguments
    args = p.parse_args()
    args.outd = args.outd.rstrip('/')
    if not os.path.exists(args.sumstats):
        sys.exit(f"Error: '{args.sumstats}' not found.")
    if not os.path.exists(args.outd):
        sys.exit(f"Error: '{args.outd}' not found.")
    if args.number_sites <= 0:
        sys.exit(f"Error: 'number-sites' ({args.number_sites}) must be non-zero positive integer.")
    # Maf and Mac cannot be specified together
    if args.min_maf is not None and args.min_mac is not None:
        sys.exit(f"Error: min_maf and min_mac are mutualy exclusive and cannot be specified together.")
    # min_num_indv and min_perc_indv cannot be specified together
    if args.min_perc_indv is not None and args.min_num_indv is not None:
        sys.exit(f"Error: min_perc_indv and min_num_indv are mutualy exclusive and cannot be specified together")
    return args

# Site from the Sumstats File
class SumstatsSite:
    def __init__(self, locus_id, locus_col, chrom, bp, popid, p, hwe, private):
        self.locid = locus_id
        self.locol = locus_col
        self.chrom = chrom
        self.bp    = bp
        self.popid = popid
        self.p     = p
        self.hwe   = hwe
        self.priv  = private
    def __str__(self):
        return f'{self.locid}\t{self.locol}\t{self.chrom}\t{self.bp}\t{self.popid}\t{self.p}\t{self.hwe}\t{self.priv}'

# Filter parameters
class FiltParams:
    def __init__(self, n_pops=1, p_inds=None, n_inds=None, hwe=False, hwe_alpha=0.05, min_maf=None, min_mac=None, max_obs_het=1.0):
        # Check inputs
        assert type(n_pops) in [int, float]
        assert n_pops >= 1, f'n_pops ({n_pops}) must be greater or equal than 1.'
        # Proportion Inds
        if p_inds is not None:
            try:
                p_inds = float(p_inds)
            except ValueError:
                sys.exit(f'Error: min_perc_indv ({p_inds}) not a floating point number.')
            assert 0.0 < p_inds <= 1.0, f'min_perc_indv ({p_inds}) must be between 0 and 1.'
        # Number Inds
        if n_inds is not None:
            try:
                n_inds = int(n_inds)
            except ValueError:
                sys.exit(f'Error: min_num_indv ({n_inds}) not an integer.')
            assert n_inds >= 1, f'n_inds ({n_inds}) must be greater or equal than 1.'
        # HWE
        assert type(hwe) is bool
        assert type(hwe_alpha) is float
        assert 0.0 < hwe_alpha <= 1.0, f'HWE alpha ({hwe_alpha}) must be between 0 and 1'
        # MAF
        if min_maf is not None:
            try:
                min_maf = float(min_maf)
            except ValueError:
                sys.exit(f'Error: min_maf ({min_maf}) not a floating point number.')
            assert 0.0 <= min_maf <= 0.5, f'Min MAF ({min_maf}) must be between 0.0 and 0.5.'
        # MAC
        if min_mac is not None:
            try:
                min_mac = int(min_mac)
            except ValueError:
                sys.exit(f'Error: min_mac ({min_mac}) not an integer.')
            assert 0 <= min_mac, f'Min MAC ({min_mac}) must be greater than 0.'
        assert type(max_obs_het) is float
        assert 0.0 <= max_obs_het <= 1.0, f'Max Observed Het ({max_obs_het}) must be between 0 and 1.'
        # Assign parameters
        self.pops  = n_pops
        self.pinds = p_inds
        self.ninds = n_inds
        self.hwe   = hwe
        self.alpha = hwe_alpha
        self.maf   = min_maf
        self.mac   = min_mac
        self.het   = max_obs_het
    def filter_minor_allele(self, p, n_indv):
        remove = False
        # Filters not applied
        if self.maf is None and self.mac is None:
            remove = False
        # Filter by MAF
        elif self.maf is not None and self.mac is None:
            if p >= 0.5:
                if 0 < 1-p < self.maf:
                    remove = True
            else:
                if 0 < p < self.maf:
                    remove = True
        # Filter by MAC
        else:
            c = p*(n_indv*2)
            if p >= 0.5:
                c = (1-p)*(n_indv*2)
            if 0 < c < self.mac:
                remove = True
        return remove
    def filter_by_indv(self, n_obs_indv, total_indv_in_pop):
        remove = False
        # Filters not applied
        if self.ninds is None and self.pinds is None:
            remove = False
        # Filter by number of individuals
        elif self.ninds is not None and self.pinds is None:
            if n_obs_indv < self.ninds:
                remove = True
        # Filter by proportion of individuals
        elif self.pinds is not None and self.ninds is None:
            p_obs_indv = n_obs_indv / total_indv_in_pop
            if p_obs_indv < self.pinds:
                remove = True
        else:
            remove = False
        return remove
    def print_to_log(self):
        log_str = 'Filtering SUMSTATS using the following filters:\n'
        log_str += f'    Min number of populations: {self.pops}\n'
        if self.pinds is not None:
            log_str += f'    Min proportion of individuals per pop: {self.pinds:.6g}\n'
        if self.ninds is not None:
            log_str += f'    Min number of individuals per pop: {self.ninds}\n'
        if self.hwe:
            log_str += f'    Removing loci out of HWE per pop (p-value cutoff {self.alpha:.6f})\n'
        if self.maf is not None:
            log_str += f'    Min minor allele frequency (MAF): {self.maf:.6g}\n'
        if self.mac is not None:
            log_str += f'    Min minor allele count (MAC): {self.mac}\n'
        if self.het < 1:
            log_str += f'    Max observed heterozygosity per pop: {self.het:.6g}\n'
        print(log_str)

#
# Print date-time
#
def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

#
# Generate the class with the filtering options
#
def parse_filter_params(prog_opt):
    filt_opts = FiltParams(
        prog_opt.min_pops,
        prog_opt.min_perc_indv,
        prog_opt.min_num_indv,
        prog_opt.hwe,
        prog_opt.hwe_alpha,
        prog_opt.min_maf,
        prog_opt.min_mac,
        prog_opt.max_obs_het)
    filt_opts.print_to_log()
    return filt_opts

#
# Parse Sumstats file
#
def parse_sumstats(sumstats_f, filt_opts):
    assert isinstance(filt_opts, FiltParams)
    # Print to log
    sum_f = sumstats_f.split('/')[-1]
    print(f'Processing {sum_f}...')
    # Output
    sites_dict = dict()
    prev_snp = None
    pops_in_site = list()
    # Determine number of samples per pop
    pop_samples = dict()
    tot_samples = 0
    # Tally
    record = 0
    loci = set()
    snps = set()
    # Parse the sumstats
    for i, line in enumerate(open(sumstats_f, 'r')):
        if line.startswith('#'):
            line = line.lstrip('#').lstrip(' ')
            if line.startswith('Locus ID'):
                # Skip the true header
                continue
            else:
                fields = line.strip('\n').split('\t')
                pop = fields[0]
                inds = fields[1].split(',')
                pop_samples[pop] = len(inds)
                tot_samples += len(inds)
                continue
        # For variant site rows
        # First, check if you got the sample IDs from the header
        assert tot_samples > 0, f'Error: No sample/population IDs were detected in the header'
        # Parse the rows
        fields = line.strip('\n').split('\t')
        locid  = int(fields[0])      # Locus ID
        locol  = int(fields[3])      # SNP column
        chrom  = fields[1]           # Chromosome ID
        bp     = int(fields[2])      # Basepair
        popid  = fields[4]           # Population ID
        n_indv = int(fields[7])      # Number of individuals 
        p_freq = float(fields[8])    # P allele freq
        ob_het = float(fields[9])    # Observed Heterozygosity
        hwe_p  = float(fields[19])   # HWE p-value
        priv   = int(fields[20])     # Private allele status
        snp    = f'{locid}_{locol}'  # SNP ID
        record += 1
        loci.add(locid)
        snps.add(snp)
        # 1. Filter by the minimum number/proportion of individuals
        pop_indv = pop_samples.get(popid, None)
        if pop_indv is None:
            sys.exit(f'Error: No samples for population \'{popid}\' seen in header.')
        # Calculate the proportion of individuals at the locus
        rm_by_inds = filt_opts.filter_by_indv(n_indv, pop_indv)
        if rm_by_inds:
            continue
        # 2. Filter by maf/mac
        rm_by_allele = filt_opts.filter_minor_allele(p_freq, n_indv)
        if rm_by_allele:
            continue
        # 3. Filter by Obs Het
        if ob_het > filt_opts.het:
            continue
        # 4. Filter by HWE
        if filt_opts.hwe:
            if hwe_p == -1.0:
                # The HWE was not used in populations
                warnings.warn('HWE P-value is \'-1.00000\'. The --hwe was not added to populations. Skipping HWE filtering.')
                filt_opts.hwe = False
                hwe_p = 1.0
            if hwe_p < filt_opts.hwe_alpha:
                continue
        # Process remaining sites
        site = SumstatsSite(locid, locol, chrom, bp, popid, p_freq, hwe_p, priv)
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
    print(f'    Read {record:,} total records from input SUMSTATS.')
    print(f'    Processed {len(loci):,} total input loci, composed of {len(snps):,} total SNPs.\n')
    print_tally(sites_dict)
    return sites_dict

#
# Print a tally of filtered loci
#
def print_tally(sites_dict):
    # Tally total kept sites
    nlocs = len(sites_dict)
    nsnps = 0
    for loc in sites_dict:
        for site in sites_dict[loc]:
            nsnps += 1
    print('After filtering, the program kept:')
    print(f'    Total loci: {nlocs:,}')
    print(f'    Total SNPs: {nsnps:,}\n')

#
# Sample the final sites
#
def sample_kept_sites(locus_dict, n_sites, outd='./', single_snp=False):
    n_loci = list(locus_dict.keys())
    export = 0
    if len(n_loci) == 0:
        return
    with open(f'{outd}/snp_whitelist.tsv', 'w') as fh:
        # Exporting single snp per-locus
        if single_snp:
            print('Exporting sites (single SNP per locus)...')
            if n_sites > len(n_loci):
                print(f'    Warning: more sites chosen for export ({n_sites:,}) than loci available post filtering ({len(n_loci):,}). Exporting {len(n_loci):,} total sites.')
                n_sites = len(n_loci)
            sampled_loci = random.sample(n_loci, n_sites)
            for locus in sorted(sampled_loci):
                cols = list(locus_dict[locus].keys())
                col = random.choice(cols)
                fh.write(f'{locus}\t{col}\n')
            print(f'    Exported {n_sites:,} total sites to the whitelist.')
        # When exporting multiple SNPs per locus:
        else:
            print('Exporting sites (more than one variant sites can be exported perlocus)...')
            # Make a new list of SNP ids, regardless of locus
            snps = list()
            for loc in sorted(locus_dict):
                cols = locus_dict[loc]
                for col in sorted(cols):
                    snp_id = (loc, col)
                    snps.append(snp_id)
            # Adjust if the user asks for more sites than available
            if n_sites > len(snps):
                print(f'    Warning: more sites chosen for export ({n_sites:,}) than sites available post filtering ({len(snps):,}). Exporting {len(snps):,} total sites.')
                n_sites = len(snps)
            sampled_snps = random.sample(snps, n_sites)
            for snp in sorted(sampled_snps):
                fh.write(f'{snp[0]}\t{snp[1]}\n')
            print(f'    Exported {len(sampled_snps):,} total sites to the whitelist.')
    fh.close()


def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    # Process input filter parameters
    filt_params = parse_filter_params(args)
    # Read sumstats and filter sites
    sites_dict = parse_sumstats(args.sumstats, filt_params)
    # Subsample the kept sites and save to file
    sample_kept_sites(sites_dict, args.number_sites, args.outd, args.write_random_snp)
    print(f'\n{PROG} finished on {now()}')

# Run Code
if __name__ == '__main__':
    main()
