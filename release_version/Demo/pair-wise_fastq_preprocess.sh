#!/bin/bash
#single known reference

input_folder="../../input"
output_folder="../../output"
input_reads_file="SRR3382456"

./main.sh -p preprocess -w unzip -r $input_folder/${input_reads_file}_1.fastq.gz -o $output_folder/pre_process
./main.sh -p preprocess -w unzip -r $input_folder/${input_reads_file}_2.fastq.gz -o $output_folder/pre_process
./main.sh -p preprocess -w merge -r $output_folder/pre_process/${input_reads_file}_1.fastq,$output_folder/pre_process/${input_reads_file}_2.fastq -o $output_folder/pre_process

# After finish preprocess step, then you can use single_known or multiple_known or species_unknown scripts to run..
# just needs to make sure when use: extract -r $input_folder/${input_reads_file}.fastq , input_folder=$output_folder/pre_process
