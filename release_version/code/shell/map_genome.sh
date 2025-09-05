#!/bin/bash

database_name=$1
qcovhsp=$2
pident=$3
blast_fna=$4
clean_fa=$5
database=$6
results=$7


echo $qcovhsp
echo $pident
# Function to check if a value is an integer
is_integer() {
    local value=$1
    [[ $value =~ ^[0-9]+$ ]]
}

# Validate qcovhsp and pident
if ! is_integer "$qcovhsp" || ! is_integer "$pident"; then
    printf "Error: qcovhsp (%s) or pident (%s) is not a valid integer\n" "$qcovhsp" "$pident" >&2
    exit 1
fi

###start: 1st step, make database
# find $blast_refprok_mapping_fna -type f -exec cat {} + > $blast_refprok_mapping_fna/combined_gcf.fna
makeblastdb -in $blast_fna -dbtype nucl -out ${database}/${database_name}

# makeblastdb -in $blast_fna -dbtype nucl -out ${database}/${database_name}
# seqtk seq -a $results/run_result/final_seq12.fastq > $results/run_result/final_seq12.fa 
blastn -db ${database}/${database_name} -query $clean_fa -num_threads 4 -task blastn -dust no -outfmt "6 delim=, qseqid sacc sstart send evalue bitscore qcovhsp pident" -max_target_seqs 5 -max_hsps 3 -out $results/blast_score.txt

# wc -l $blast_score/blast_score_1db_${bacteria}.txt
# awk -F, '{print $1}' $blast_score/blast_score_1db_${bacteria}.txt | sort | uniq -c | wc -l
# awk -F, '($7 >= 80 && $8 >= 90) {print $1}' $blast_score/blast_score_1db_${bacteria}.txt | sort | uniq -c | wc -l
awk -F, -v pq=$qcovhsp -v pp=$pident '($7 >= pq && $8 >= pp) {print $0}' $results/blast_score.txt  > $results/blast_score_filter.txt
# awk -F, '($7 >= $qcovhsp && $8 >= $pident) {print $0}' $results/blast_score.txt  > $results/blast_score_filter.txt

number=$(awk -F, '($7 >= 80 && $8 >= 90) {print $1}' $results/blast_score_filter.txt | sort | uniq -c | wc -l)
number_2=$(awk -F, '($7 >= 90 && $8 >= 90) {print $1}' $results/blast_score_filter.txt | sort | uniq -c | wc -l)
number_3=$(awk -F, '($7 >= 100 && $8 >= 90) {print $1}' $results/blast_score_filter.txt | sort | uniq -c | wc -l)
number_4=$(awk -F, '($7 >= 100 && $8 >= 100) {print $1}' $results/blast_score_filter.txt | sort | uniq -c | wc -l)
echo "Based on the paramter qcovhsp:$qcovhsp and pident:$pident you set, we save blast_score_fileter.txt based on this creteria." >> $results/genome_mapping_analysis.csv
total_number=$(wc -l < $clean_fa | awk '{print $1/2}')
echo "Except this, we can give you the standard mapping results overview as follows:" >> $results/genome_mapping_analysis.csv
echo "mapping top species database number(qcovhsp=80,pident=90): $number" >> $results/genome_mapping_analysis.csv
echo "mapping top species database number(qcovhsp=90,pident=90): $number_2" >> $results/genome_mapping_analysis.csv
echo "mapping top species database number(qcovhsp=100,pident=90): $number_3" >> $results/genome_mapping_analysis.csv
echo "mapping top species database number(qcovhsp=100,pident=100): $number_4" >> $results/genome_mapping_analysis.csv
echo "total sequence number(after clean): $total_number" >> $results/genome_mapping_analysis.csv
echo "$number $total_number" | awk '{print "mapping percentage(qcovhsp=80,pident=90):", ($1/$2)*100 "%"}' >> $results/genome_mapping_analysis.csv
echo "$number_2 $total_number" | awk '{print "mapping percentage(qcovhsp=90,pident=90):", ($1/$2)*100 "%"}' >> $results/genome_mapping_analysis.csv
echo "$number_3 $total_number" | awk '{print "mapping percentage(qcovhsp=100,pident=90):", ($1/$2)*100 "%"}' >> $results/genome_mapping_analysis.csv
echo "$number_4 $total_number" | awk '{print "mapping percentage(qcovhsp=100,pident=100):", ($1/$2)*100 "%"}' >> $results/genome_mapping_analysis.csv

# cp $blast_score/blast_score_filter.txt $blast_refprok_mapping_results/blast_score_filter_1db_${1}.txt