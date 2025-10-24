#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Global variables
CODE_PATH_AFTER_CLEAN=$1
FINAL_FASTA=$2
MAPPING_MIRNA_FILE=$3
REFERENCE_GENOME_FILE=$4
MAPPING_GENOME_FILE=$5
DUPLICATES_INFORMATION=$6
MIRDEEP_FOLDER=$7
LINEARFOLD_RESULT=$8
OUTPUT_FOLDER=$9
echo "a$OUTPUT_FOLDER"
# Function to produce overlap analysis
produce_final_form() {
    local final_fasta=$1
    local mapping_mirna_file=$2
    local reference_genome_file=$3
    local mapping_genome_file=$4
    local duplicates_information=$5
    local mirdeep_folder=$6
    local linearfold_result=$7
    local output_folder=$8
    
    python3 "${CODE_PATH_AFTER_CLEAN}/produce_final_form.py" \
        -input_final_fasta "${final_fasta}" \
        -input_genome_mapping_results "${mapping_genome_file}" \
        -input_reference_genome "${reference_genome_file}" \
        -duplicates_information "${duplicates_information}" \
        -input_mirna_mapping_results "${mapping_mirna_file}" \
        -input_mirdeep2_folder "${mirdeep_folder}" \
        -input_linearfold_results "${linearfold_result}" \
        -output_folder "${output_folder}" 
}

produce_final_form "${FINAL_FASTA}" "${MAPPING_MIRNA_FILE}" "${REFERENCE_GENOME_FILE}" "${MAPPING_GENOME_FILE}" "${DUPLICATES_INFORMATION}" "${MIRDEEP_FOLDER}" "${LINEARFOLD_RESULT}" "${OUTPUT_FOLDER}"