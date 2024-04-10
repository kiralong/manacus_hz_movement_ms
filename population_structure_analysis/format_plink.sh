#!/bin/bash

working_dir=/path/to/working/directory/admixture
old_ped=$working_dir/populations.plink.ped
new_ped=$working_dir/populations.corrected.plink.ped
old_map=$working_dir/populations.plink.map
new_map=$working_dir/populations.corrected.plink.map

# This series of sed commands renames the Manacus sampling site name to be numeric only
cat $old_ped | grep -v "^#" | \
    sed 's/^pr_090/090/' | \
    sed 's/^ru_080/080/' | \
    sed 's/^cg_100/100/' | \
    sed 's/^ss_020/020/' | \
    sed 's/^qp_060/060/' | \
    sed 's/^ro_050/050/' | \
    sed 's/^fc_040/040/' | \
    sed 's/^mr_095/095/' | \
    sed 's/^so_030/030/' > $new_ped

cat $old_map | grep -v "^#" > $new_map

# Plink command to convert .ped and .map files to bed file
plink --file populations.corrected.plink --make-bed --recode12 --allow-extra-chr 0 --out populations.corrected.plink
