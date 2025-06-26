#!/bin/bash
# Exit immediately if any command fails
source ~/.bashrc
conda activate cleanfasq_env
set -e

# Function to delete file if it exists
delete_file_if_exists() {
    local file_path="$1"
    if [[ -f "$file_path" ]]; then
        rm -f "$file_path"
    fi
}

update_combined_fna_format() {
    local input_file="$1"
    local temp_file="${input_file}.tmp"
    local header=">chr_default_combined:1-99999"

    if [[ -f "$input_file" ]]; then
        grep -v '^>' "$input_file" > "$temp_file"
        printf "%s\n" "$header" | cat - "$temp_file" > "$input_file"
        rm -f "$temp_file"
    else
        printf "Error: File %s not found\n" "$input_file" >&2
        return 1
    fi
}


database_name=$1
code_path=$2
clean_fasta_file=$3
fna_file=$4
results=$5


echo "start"
echo $(date)
# ## nead filter original reads only match species part
# python3 ${code_path}/prepare_mirdeep_reads_file.py -input_reads $clean_fasta_file -csv_file $mapping_filter_score_species_file -output_file $results/middle_results/reads.fa
cp $clean_fasta_file $results/middle_results/reads.fa
echo $(date)

## start mirdeep2 preparation
# delete_file=$results/middle_results/reads_collapsed_vs_genome.arf
delete_file_if_exists "$results/middle_results/reads_collapsed_vs_genome.arf"
cp $fna_file $results/middle_results/fna_file.fa
update_combined_fna_format "$results/middle_results/fna_file.fa"

cd $results
## mkdir -p $results/middle_results/database

## bowtie to generate database
module load bowtie/1.3
bowtie-build $results/middle_results/fna_file.fa $results/middle_results/database/${database_name}

## start mirdeep2
set +e
collapse_reads_md.pl $results/middle_results/reads.fa seq > $results/middle_results/reads_collapsed.fa
mapper.pl $results/middle_results/reads_collapsed.fa -c -p $results/middle_results/database/${database_name} -t $results/middle_results/reads_collapsed_vs_genome.arf
miRDeep2.pl $results/middle_results/reads_collapsed.fa $results/middle_results/fna_file.fa $results/middle_results/reads_collapsed_vs_genome.arf none none none 2>$results/report.log
set -e
## end