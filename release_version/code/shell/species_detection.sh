#!/bin/bash
# Exit immediately if any command fails
set -e

refprok_database=$1
code_path=$2
top_num=$3
clean_fasta=$4
results=$5
temp_folder=$results/middle_results
blast_refprok_mapping_fna=$temp_folder/ncbi_dataset_fna/ncbi_dataset/data

echo "$refprok_database"
echo "$code_path"
echo "$clean_fasta"
echo "$results"

module load blast/2.14
module load biodata/1.0

##start, 1st step
# blastn -query $clean_fasta -db ${refprok_database}/ref_prok_rep_genomes -out $results/middle_results/blastn_refprok_evalue10.txt -num_threads 4 -evalue 10 -outfmt "6 qseqid sacc sstart send evalue bitscore qcovhsp pident sseqid staxids sscinames sblastnames" -max_target_seqs 5 -max_hsps 3 -task "blastn"

# wc -l $results/middle_results/blastn_refprok_evalue10.txt
# awk -F'\t' '{print $1}' $results/middle_results/blastn_refprok_evalue10.txt | sort -n | uniq -c | wc -l

###2nd step
num=$(wc -l $clean_fasta |awk '{print $1/2}')
echo $num
python3 ${code_path}/blast_species_classification.py -input_folder $results/middle_results -total_sequence_number $num  -output_folder $results

###3rd step
python3 ${code_path}/union_top_species_classification.py -input_folder $results/middle_results -top_count $top_num -output_folder $results

###4th step: this step need user add all not found gcf_numbers to csv
##hand made step

###optional 5th step
gcf_numbers=$(awk -F'\t' 'NR > 1 && $2 {print $2}' $results/mapping.csv | tr '\n' ' ')
echo $gcf_numbers

curl -o $temp_folder/datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
chmod 777 $temp_folder/datasets
${temp_folder}/datasets download genome accession ${gcf_numbers} --include genome --filename $temp_folder/ncbi_dataset_fna.zip
# ${temp_folder}/datasets download genome accession ${gcf_numbers} --include genome --format gtf --filename $temp_folder/ncbi_dataset_gtf.zip
unzip -o $temp_folder/ncbi_dataset_fna.zip -d $temp_folder/ncbi_dataset_fna
find $temp_folder/ncbi_dataset_fna/ncbi_dataset/data -type f -name "*.fna" -exec cat {} + > $results/combined_gcf.fna

