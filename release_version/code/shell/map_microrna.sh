#!/bin/bash

database_name=$1
code_path=$2
fa_file=$3
clean_fasta_file=$4
database=$5
results=$6
# Be sure to use a directory under $WORK for your job

makeblastdb -in $fa_file -dbtype nucl -out ${database}/${database_name}
blastn -query $clean_fasta_file -db ${database}/${database_name} -out $results/blastn_hairpin_rna.txt -word_size 13 -gapopen 5 -gapextend 2 -num_threads 4 -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 10 -evalue 100 -strand "plus" -task blastn
###step2: calculate percentage
num=$(wc -l $clean_fasta_file |awk '{print $1/2}')
python3 $code_path/calculate_miRNA_mapping_percentage.py -total_seq_num $num -input_path $results/blastn_hairpin_rna.txt -output_folder $results

python3 $code_path/extract_hairpin_mapping_sequences.py -input_fasta $clean_fasta_file -input_mapping $results/blastn_hairpin_rna.txt -output_path $results/blastn_hairpin_sequences.csv
