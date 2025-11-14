#!/bin/bash
#single known reference

input_folder="../../input"
output_folder="../../output"
fna_folder="../../database/fna_folder/dataset1"
gtf_folder="../../database/gtf_folder/dataset1"
hairpin_folder="../../database/hairpin_folder"
input_reads_file="SRR18745680"
tag="*AGATCGGAAGAGCACA*"
length="12"
tolerance_1="2"
tolerance_2="4"
umi_parameter="0"
qcovhsp="100"
pident="100"
bacteria="Pseudomonas_aeruginosa_PAO1"
gcf_number="GCF_000006765.1_ASM676v1_genomic"
umi_file="none"
shortest_sequence_length="18"
longest_sequence_legnth="40"

./main.sh -p extract -r $input_folder/${input_reads_file}.fastq -o $output_folder/extract -l $length -F $tag  --t1 $tolerance_1 --t2 $tolerance_2 --umi $umi
./main.sh -p map_genome -f $fna_folder${gcf_number}.fna -c $output_folder/extract/final_seq_${length}.fa -d $output_folder/map_genome/${bacteria} -n ${bacteria}  --pq $qcovhsp --pp $pident -o $output_folder/map_genome
./main.sh -p map_mirna -f $hairpin_folder/hairpin.fa -c $output_folder/extract/final_seq_${length}.fa -d $output_folder/map_mirna/hairpin_database -n hairpin  -o $output_folder/map_mirna
./main.sh -p quantify_genome -g $gtf_folder/${gcf_number}.gtf -m $output_folder/map_genome/blast_score_filter.txt -u $umi_file -o $output_folder/quantify_genome
./main.sh -p quantify_mirna  -m $output_folder/map_mirna/blastn_hairpin_rna.txt -u $umi_file -o $output_folder/quantify_mirna
./main.sh -p integrate -c $output_folder/extract/final_seq_${length}.fastq -s $output_folder/map_genome -m $output_folder/map_mirna -u $umi_file --sl $shortest_sequence_length --ll $longest_sequence_legnth -o $output_folder/integrate
./main.sh -p predict -w mirdeep2 -c $output_folder/integrate/prediction_input.fasta  -f $fna_folder${gcf_number}.fna -o $output_folder/predict_mirdeep -n database
./main.sh -p predict -w linearfold  -c $output_folder/integrate/prediction_input.fasta -m $output_folder/integrate/blast_score_prediction_filter.txt   -f $fna_folder${gcf_number}.fna  -o $output_folder/predict_linearfold
./main.sh -p produce_final_form -c $output_folder/integrate/prediction_input.fasta -m $output_folder/map_mirna/blastn_hairpin_sequences.csv --inf $output_folder/integrate/redundant_sequences_information.csv --mr $output_folder/predict_mirdeep --lf $output_folder/predict_linearfold/hairpin_information.csv  -o $output_folder/produce_final_form

