#!/bin/bash
# Exit immediately if any command fails
set -e

code_path=$1
filter_score_file=$2
umi_file=$3
results=$4


cp ${filter_score_file} ${results}/middle_results/blastn_hairpinrna.txt
###miRNA quantification, no need gtf data
python3 ${code_path}/blast_process_main.py find_unique_count_mirna -input_sorted_score_file $results/middle_results/blastn_hairpinrna.txt -input_umi_file $umi_file -output_unique ${results}/output_unique_biotype_mirna.csv
