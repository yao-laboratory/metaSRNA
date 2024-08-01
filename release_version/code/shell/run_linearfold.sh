#!/bin/bash

source ~/.bashrc
conda activate cleanfasq_env
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
cat ${middle_results}/predict_use.fasta | $linearfold_folder/./linearfold --fasta --Vienna --verbose | grep 'index\|Hairpin\|\.\.' > ${middle_results}/predict_final
echo $(date)

# ### start linear fold result analysis
# export PARAM1=$results/predict_result/${sample_id}_my_predict
# export PARAM2=$results/blast_result/${sample_id}_seq12_unique_sequences.fasta
# export PARAM3=$results/blast_result/blastn_filter_${sample_id}_shorter_after_filter.csv
# export PARAM4=$results/predict_result/${sample_id}_hairpin_information.csv
echo $(date)
python3 $final_linear_fold_result ${middle_results}/predict_final $clean_fasta_file $mapping_filter_score_species_file ${results}/hairpin_information.csv
echo $(date)