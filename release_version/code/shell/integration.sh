#!/bin/bash

module load blast/2.14
module load biodata/1.0

code_path_for_clean=$1
code_path_after_clean=$2
input_fastq=$3
mapping_species=$4
mapping_mirna=$5
results=$6

python3 $code_path/filter_blast_results_length.py -input_file $mapping_species -start 18 -end 40 -output_folder $results
python3 ${code_path_for_clean}/process_sequence_length.py process_length -input $input_fastq -filter_min_length 18 -filter_max_length 40 -output_folder $results/middle_results
seqtk seq -a $results/middle_results/final_seq_18_to_40.fastq  > $results/middle_results/final_seq_18_to_40.fasta

## filter species mapping file with new reads length requirements
python3 ${code_path_after_clean}/filter_blast_result_length.py -input_mapping_file $mapping_species/blast_score_filter.txt -input_filtered_clean_fasta $results/middle_results/final_seq_18_to_40.fasta -output_file ${results}/blast_score_final_filter.txt
### produce mapping species file and mapping mirna file venn diagram results
python3 ${code_path_after_clean}/mirna_and_species_blast_overlap.py -species_mapping_file ${results}/blast_score_final_filter.txt -hairpin_folder $mapping_mirna -overlap_folder $results

# ## nead filter filtering reads only match species part, this is an input for predict step
python3 ${code_path_after_clean}/prepare_reads_file.py -input_reads $results/middle_results/final_seq_18_to_40.fasta -csv_file $mapping_species/blast_score_filter.txt -output_file $results/final_species.fasta


