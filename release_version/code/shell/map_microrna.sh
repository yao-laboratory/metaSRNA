#!/bin/bash

module load blast/2.14
module load biodata/1.0
module load allinea
module load seqtk

fa_file=$1
clean_fasta_file=$2
database=$3
database_name=$4
code_path=$5
results=$6
# Be sure to use a directory under $WORK for your job

makeblastdb -in $fa_file -dbtype nucl -out ${database}/${database_name}
blastn -query $clean_fasta_file -db ${database}/${database_name} -out $results/blastn_hairpin_rna.txt -word_size 13 -gapopen 5 -gapextend 2 -num_threads 4 -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 10 -evalue 100 -strand "plus" -task blastn
###step2: calculate percentage
num=$(wc -l $clean_fasta_file |awk '{print $1/2}')
python3 $code_path/calculate_miRNA_mapping_percentage.py -total_seq_num $num -input_path $results/blastn_hairpin_rna.txt -output_folder $results
