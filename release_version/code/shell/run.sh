#!/bin/bash
# Exit immediately if any command fails
set -e

input=$1
results=$2
code_path=$3
echo "$results"
echo "$code_path"

# clean_command_py=${code_path}/clean_command.py
# process_length_py=${code_path}/process_sequence_length.py


# gzip -d -c $input/$1.fastq.gz > $results/$1.fastq
python3 ${code_path}/clean_command.py clean_fastq -input $input -output_filename $results/middle_results/output -fa_format *TGGAATTCTCGGGTGCCAAGGAACTCCA*  -fault_tolerance 2 -tail_incomplete_tolerance 4
cat $results/middle_results/output_1_step1.fastq $results/middle_results/output_1_step2.fastq > $results/middle_results/final_seq.fastq 
cat $results/middle_results/output_2_step1.fastq $results/middle_results/output_2_step2.fastq > $results/final_umi.fastq 
python3 ${code_path}/process_sequence_length.py process_length -input $results/middle_results/final_seq.fastq  -output_folder $results -filter_length 12

seqtk seq -a $results/final_seq12.fastq  > $results/final_seq12.fasta
seqtk seq -a $results/final_seq12.fastq  > $results/final_seq12.fa
# wc -l $results/$1.fastq | awk '{print $1/4}'
# wc -l $results/run_result/final.fastq | awk '{print $1/4}'
# wc -l $results/run_result/final_seq12.fastq | awk '{print $1/4}'