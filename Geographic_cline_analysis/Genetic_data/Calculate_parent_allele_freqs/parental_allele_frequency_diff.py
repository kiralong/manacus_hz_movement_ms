populations_sumstats_file_path = "data/long_populations.sumstats.tsv"
snp_cline_parameters_file_path = "data/long_snp_cline_parameters.filtered.tsv"
out_snp_cline_parameters_file_path = "out.snp_cline_parameters_with_parental_allele_freq_diff.tsv"

# Parental Populations of interest - the two non-hybrid parental populations. 
# NOTE: Script only supports two parental populations.
parental_populations = ('cg_100', 'ss_020')


filtered_population_sumstats = []
with open(populations_sumstats_file_path, "r") as populations_sumstats_file:
    lines=populations_sumstats_file.readlines()
    
    for i in range(0,len(lines)):
        line = lines[i]

        # Get rid of headers -  If comment line, skip
        if line.startswith('#'):
            continue

        [locus_id, chr, bp, col, pop_id, p_nuc, q_nuc, n, p, *_ ] = line.split()
        if pop_id in parental_populations:
            # Create the snp_id "<Locus ID>_<Col>"
            filtered_population_sumstats.append((f"{locus_id}_{col}", p))

snp_parental_diffs = []
for i in range(0, len(filtered_population_sumstats), 2):
    p1 = filtered_population_sumstats[i]
    p2 = filtered_population_sumstats[i+1]
    diff = abs(float(p1[-1]) - float(p2[-1]))
    snp_parental_diffs.append((p1[0], diff))

with open(snp_cline_parameters_file_path, "r") as snp_cline_parameters_file:
    snp_cline_parameters_lines = snp_cline_parameters_file.readlines()

# FASTER ALGORITHM (which is still slow)
# Since the snp_clines params are a subset of the sumstats
# and thus are considerably smaller, 
# iterate through the snp_cline_parameters, scanning through the sumstats
# This lower-bounds the algorithm to the size of cline params which is smaller.

print(f"Processing snp_cline_parameters_lines {len(snp_cline_parameters_lines)} x snp_parental_diffs {len(snp_parental_diffs)}", )
for i in range(1, len(snp_cline_parameters_lines)):
    if i % 1000 == 0:
        print("PROGRESS snp_cline_parameters_lines", f"i={i}/{len(snp_cline_parameters_lines)}")

    snp_cline_parameters_line = snp_cline_parameters_lines[i]
    [_, snp_id, *_ ] = snp_cline_parameters_line.split()

    found = False
    for j in range(0, len(snp_parental_diffs)):
        snp_parental_diff = snp_parental_diffs[j]        

        if snp_parental_diff[0] == snp_id:
            snp_cline_parameters_line = snp_cline_parameters_line.rstrip() + f"\t{snp_parental_diff[1]}\n"
            snp_cline_parameters_lines[i] = snp_cline_parameters_line
            found = True
            continue
    
    if not found:
        # This should technically never happen... if it does, there's a problem.
        print("Found an NA", f"i={i} j={j}", snp_cline_parameters_line)
        snp_cline_parameters_line = snp_cline_parameters_line.rstrip() + f"\tNA\n"
        snp_cline_parameters_lines[i] = snp_cline_parameters_line

# One last thing: Add a column header on the first line
snp_cline_parameters_lines[0] = snp_cline_parameters_lines[0].rstrip() + f"\tparent_allele_freq_diff\n"

# Output the cline parameters file
with open(out_snp_cline_parameters_file_path, "w") as outfile:
    outfile.writelines(snp_cline_parameters_lines)
