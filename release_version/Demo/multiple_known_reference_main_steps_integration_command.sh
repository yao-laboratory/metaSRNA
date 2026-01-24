#!/bin/bash
# multiple known reference
# needs combine difference references fna to combined.fna and  combine difference references gtf to combined.gtf

output_folder="../../output"
fna_folder="../../database/fna_folder/dataset11"
gtf_folder="../../database/gtf_folder/dataset11"
hairpin_folder="../../database/hairpin_folder"
input_reads_file="SRR684065"
tag="*TCGTAT*"
length="12"
tolerance_1="2"
tolerance_2="4"
umi_parameter="0"
qcovhsp="100"
pident="100"
bacteria="Escherichia_coli_group"
gcf_number="combined"
umi_file="none"
shortest_sequence_length="18"
longest_sequence_legnth="40"

./main.sh -p all -r $input_folder/${input_reads_file}.fastq -l $length -F ${tag}  -n ${bacteria} --t1 $tolerance_1 --t2 $tolerance_2 --umi $umi_parameter --pq $qcovhsp --pp $pident --sl $shortest_sequence_length --ll $longest_sequence_legnth --fna $fna_folder/${gcf_number}.fna  --gtf $gtf_folder/${gcf_number}.gtf --hairpin $hairpin_folder/hairpin.fa -o $output_folder/all_steps