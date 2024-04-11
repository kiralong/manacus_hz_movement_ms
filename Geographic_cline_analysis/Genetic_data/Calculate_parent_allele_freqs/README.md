# Allele Frequency Difference between Parents

Python script to filter for allele frequency difference between parents to certain populations.

The goal is to get the difference in parental allele frequencies and add them as a column to the cline parameters table. This will be an additional column in the snp cline parameters table.

The data in the column can then be used elsewhere, e.g. for a histogram in R.

## Installation

Install the tools described in [`.tool-versions`](./.tool-versions)

- It is recommended you install and use [`asdf`](https://asdf-vm.com/) which works directly with the `.tool-versions` file. However, it is possible to manually install all of these tools with the right versions and simply make them available on your `PATH`.

## Usage

Set the values in the variables at the start of `parental_allele_frequency.py`.

Run `parental_allele_frequency.py`:

```sh
$ python3 parental_allele_frequency.py
```

## Original Instructions

Get the difference in parental allele frequencies and add them as a column to the cline parameters table that will be graphable and sortable in R.
>
Input files:
populations.sumstats.tsv
snp_cline_parameters.filtered.tsv

Ultimately we need to add a column to the snp_cline_parameters.filtered.tsv file that has the difference in allele frequencies between the two parental populations: cg_100 and ss_020.

To calculate the allele frequency difference, you’ll need to use the 9th column labeled “P” in the sumstats.tsv file. Note the sumstats file also has a header that is denoted by # at the beginning of each line. Each snp will have 9 rows and thus 9 “P” values. We only care about the first and last P values. We want the difference between the P value for the row of cg_100 and ss_020. However, there won’t always be 9 rows there are sometimes 8 so we need to pull the rows for only cg_100 and ss_020. We also need to make sure that the difference is absolute valued so the frequency difference is always positive.

Next you need to make the snp_id from the sumstats file match what is in the cline parameters file. Each cline in the cline parameters file has a unique snp in the form of “Locus ID_Col” which was originally made from the sumstats file. Column 1 in the sumstat file is Locus ID and column 4 is Col.

Note that there are a bunch of snps in the sumstats file that will NOT be in the cline parameters file. That’s ok, but make sure the python script moves on if it doesn’t find a snp ID in the cline parameters file.

So the python script will:
Calculate the allele frequency difference between the parental pops cg_100 and ss_020
Make the unique snp ID in the sumstats file
Match the snp ID in both files
Add a column to cline parameters file with the allele frequency absolute value difference

## Authors
Kira M. Long
George W. Pantazes
