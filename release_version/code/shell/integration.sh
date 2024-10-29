#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Load required modules
module load blast/2.14
module load biodata/1.0

# Global variables
CODE_PATH_FOR_CLEAN=$1
CODE_PATH_AFTER_CLEAN=$2
INPUT_FASTQ=$3
MAPPING_SPECIES=$4
MAPPING_MIRNA=$5
UMI_FILE=$6
RESULTS=$7



# Function to process sequence length
filter_sequence_length() {
    local input_fastq=$1
    local output_folder=$2
    
    python3 "${CODE_PATH_FOR_CLEAN}/filter_sequence_length.py" process_length \
        -input "${input_fastq}" \
        -filter_min_length 18 \
        -filter_max_length 40 \
        -output_folder "${output_folder}"
}

# Function to convert fastq to fasta
convert_fastq_to_fasta() {
    local input_fastq=$1
    local output_fasta=$2

    seqtk seq -a "${input_fastq}" > "${output_fasta}"
}

# Function to prepare miRDeep2 input
prepare_prediction_input() {
    local input_fasta=$1
    local output_folder=$2
    local mapping_species=$3
    local umi_file=$4
    
    python3 "${CODE_PATH_AFTER_CLEAN}/prepare_reads_file.py" \
        -input_reads "${input_fasta}" \
        -csv_file "${mapping_species}/blast_score_filter.txt" \
        -output_file "${output_folder}/prediction_temp_input.fasta"

    # Filter redundant sequences and prepare prediction input
    python3 "${CODE_PATH_AFTER_CLEAN}/filtering_redundant_sequences.py" \
        -seq_reads "${output_folder}/prediction_temp_input.fasta" \
        -umi_reads "${umi_file}" \
        -csv_file "${mapping_species}/blast_score_filter.txt" \
        -output_fasta "${output_folder}/prediction_input.fasta" \
        -output_csv "${output_folder}/redundant_sequences_information.csv"

    # Filter species mapping file with new reads
    python3 "${CODE_PATH_AFTER_CLEAN}/filter_mapping_result.py" \
        -input_mapping_file "${mapping_species}/blast_score_filter.txt" \
        -input_filtered_clean_fasta "${output_folder}/prediction_input.fasta" \
        -output_file "${output_folder}/blast_score_prediction_filter.txt"
}


# Function to produce overlap analysis
produce_overlap_analysis() {
    local species_mapping_file=$1
    local hairpin_folder=$2
    local overlap_file=$3
    
    python3 "${CODE_PATH_AFTER_CLEAN}/mirna_and_species_blast_overlap.py" \
        -species_mapping_file "${species_mapping_file}" \
        -hairpin_folder "${hairpin_folder}" \
        -overlap_file "${overlap_file}"
}

main() {
    local middle_results="${RESULTS}/middle_results"
    # local umi_file="path_to_umi_file"  # Define the path to the umi_file

    # Process sequence length and convert to fasta
    filter_sequence_length "${INPUT_FASTQ}" "${middle_results}"
    convert_fastq_to_fasta "${middle_results}/final_seq_18_to_40.fastq" "${middle_results}/final_seq_18_to_40.fasta"

    # Prepare prediction (miRDeep2 and linearfold) input
    prepare_prediction_input "${middle_results}/final_seq_18_to_40.fasta" "${RESULTS}" "$MAPPING_SPECIES" "$UMI_FILE" 

    # Produce overlap analysis
    produce_overlap_analysis "${MAPPING_SPECIES}/blast_score_filter.txt" "${MAPPING_MIRNA}" "${RESULTS}/mirna_and_top_species_analysis.csv"
    ###still have problems, mirna mapping sequence part not unique
    # produce_overlap_analysis "${RESULTS}/blast_score_linearfold_filter.txt" "${MAPPING_MIRNA}" "${RESULTS}/mirna_and_top_species_analysis_unique_sequences.csv"
}

main "$@"
