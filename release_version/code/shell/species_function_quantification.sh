#!/bin/bash
# Exit immediately if any command fails
set -e

code_path=$1
gtf_file=$2
filter_score_file=$3
umi_file=$4
results=$5

# grep -v '^#' ${gtf_folder}/GCF_000238075.gtf > ${gtf_folder}/refine.gtf
python3 ${code_path}/blast_process_main.py add_gene -input_gcf $gtf_file -input_score $filter_score_file -output_csv ${results}/middle_results/blast_score_filter_add_gene.csv

python3 ${code_path}/blast_process_main.py find_unique_count -input_sorted_score_file ${results}/middle_results/blast_score_filter_add_gene.csv -input_umi_file $umi_file -input_keys "gene_biotype" -output_unique ${results}/output_unique_biotype.csv

python3 ${code_path}/blast_process_main.py find_unique_count -input_sorted_score_file ${results}/middle_results/blast_score_filter_add_gene.csv -input_umi_file $umi_file -input_keys "gene gene_biotype" -output_unique ${results}/output_unique_gene_biotype.csv

python3 ${code_path}/blast_process_main.py find_unique_count -input_sorted_score_file ${results}/middle_results/blast_score_filter_add_gene.csv -input_umi_file $umi_file -input_keys "gene" -output_unique ${results}/output_unique_gene.csv

