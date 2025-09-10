#!/bin/bash

# Exit immediately if any command fails
set -e

code_path=$1
clean_fasta_file=$2
mapping_filter_score_species_file=$3
fna_file=$4
results=$5
middle_results=$results/middle_results

echo $(date)
python3 $code_path/blast_process_main.py extract_sequences_from_genome -input_blast_result $mapping_filter_score_species_file -input_genome $fna_file  -output_fasta ${middle_results}/predict_use.fasta
echo $(date)

echo $(date)
cat ${middle_results}/predict_use.fasta | linearfold --fasta --Vienna --verbose | grep 'index\|Hairpin\|\.\.' > ${middle_results}/predict_final
echo $(date)

# ### start linear fold result analysis
echo $(date)
python3 $code_path/final_linear_fold_result.py -input_linearfold_result ${middle_results}/predict_final -input_reads $clean_fasta_file -csv_file $mapping_filter_score_species_file -output_file ${results}/hairpin_information.csv
echo $(date)