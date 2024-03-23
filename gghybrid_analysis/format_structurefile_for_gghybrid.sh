#!/bin/bash

structure_file=/projects/aces/kira2/stacks-2.60_run/testing_5each_pop_2022_feb_16/gghybrid/populations.structure
output_file=populations.p9.r80.sSNP.structure.csv

cat $structure_file | grep -v "^#" | tr '\t' ',' | sed -E 's/,0/,NA/g' |  sed -E 's/,,/INDLABEL,POPID,/' | \
	sed 's/,pr_090/,090PR/' | \
	sed 's/,ru_080/,080RU/' | \
	sed 's/,cg_100/,100CG/' | \
	sed 's/,ss_020/,020SS/' | \
	sed 's/,qp_060/,060QP/' | \
	sed 's/,ro_050/,050RO/' | \
	sed 's/,fc_040/,040FC/' | \
	sed 's/,mr_095/,095MR/' | \
	sed 's/,so_030/,030SO/' > $output_file
