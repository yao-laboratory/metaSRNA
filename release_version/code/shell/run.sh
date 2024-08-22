#!/bin/bash
# Exit immediately if any command fails
set -e

code_path=$1
num=$2
format=$3
fault_tl=$4
incomplete_tl=$5
input=$6
results=$7

echo "$format"
echo "$fault_tl"
echo "$incomplete_tl"
echo "$input"
echo "$results"
echo "$code_path"

# clean_command_py=${code_path}/clean_command.py
# process_length_py=${code_path}/process_sequence_length.py


# gzip -d -c $input/$1.fastq.gz > $results/$1.fastq

python3 ${code_path}/clean_command.py clean_fastq -input $input -output_filename $results/middle_results/output -fa_format $format  -fault_tolerance $fault_tl -tail_incomplete_tolerance $incomplete_tl
cat $results/middle_results/output_1_step1.fastq $results/middle_results/output_1_step2.fastq > $results/middle_results/final_seq.fastq 
cat $results/middle_results/output_2_step1.fastq $results/middle_results/output_2_step2.fastq > $results/final_umi.fastq 
# need to add: ${code_path}/clean_command.py De-duplication
python3 ${code_path}/process_sequence_length.py process_length -input $results/middle_results/final_seq.fastq  -output_folder $results -filter_min_length $num

seqtk seq -a $results/final_seq_${num}.fastq  > $results/final_seq_${num}.fasta
seqtk seq -a $results/final_seq_${num}.fastq  > $results/final_seq_${num}.fa
seqtk seq -a $results/final_umi.fastq   > $results/final_umi.fasta

# wc -l $results/$1.fastq | awk '{print $1/4}'
# wc -l $results/run_result/final.fastq | awk '{print $1/4}'
# wc -l $results/run_result/final_seq12.fastq | awk '{print $1/4}'