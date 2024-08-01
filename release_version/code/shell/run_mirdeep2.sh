#!/bin/bash
# Exit immediately if any command fails
set -e

datbase_name=$1
code_path=$2
raw_data=$3
clean_fasta_file=$4
mapping_filter_score_species_file=$5
fna_file=$6
results=$7

### need cd to the folder which will let all the running file to this folder
cd $clean_fastaq/shared/code/after_cleaning
set -e
# export PARAM1=$results/blast_result/${sample_id}_seq12_unique_sequences.fasta
# export PARAM2=$results/blast_result/blastn_filter_${sample_id}_shorter_after_filter.csv
# export PARAM3=$bowtie_mirdeep_result/reads.fa
echo "start"
echo $(date)
python3 ${code_path}/prepare_mirdeep_reads_file.py $clean_fasta_file $mapping_filter_score_species_file $$raw_data
echo $(date)

###start
cp $fna_file ${middle_results}/fna_file.fa
cp $raw_data ${middle_results}/reads.fa
mkdir -p ${middle_results}/database

cd $bowtie_mirdeep_result
module load bowtie/1.3
module load miRDeep/2.0
bowtie-build ${middle_results}/fna_file.fa ${middle_results}/database/${database_name}
collapse_reads_md.pl ${middle_results}/reads.fa mmu > $middle_results/reads_collapsed.fa
mapper.pl $middle_results/reads_collapsed.fa -c -p ${middle_results}/database/${database_name} -t $middle_results/reads_collapsed_vs_genome.arf
miRDeep2.pl $middle_results/reads_collapsed.fa ${middle_results}/fna_file.fa $middle_results/reads_collapsed_vs_genome.arf none none none 2>$results/report.log
###end