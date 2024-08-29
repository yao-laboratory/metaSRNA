#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Global variables
CODE_PATH_AFTER_CLEAN=$1
FINAL_FASTA=$2
MAPPING_MIRNA_FILE=$3
DUPLICATES_INFORMATION=$4
MIRDEEP_FOLDER=$5
LINEARFOLD_RESULT=$6
OUTPUT_FOLDER=$7
echo "a$OUTPUT_FOLDER"
# Function to produce overlap analysis
produce_final_form() {
    local final_fasta=$1
    local mapping_mirna_file=$2
    local duplicates_information=$3
    local mirdeep_folder=$4
    local linearfold_result=$5
    local output_folder=$6
    
    python3 "${CODE_PATH_AFTER_CLEAN}/produce_final_form.py" \
        -input_final_fasta "${final_fasta}" \
        -duplicates_information "${duplicates_information}" \
        -input_mirna_mapping_results "${mapping_mirna_file}" \
        -input_mirdeep2_folder "${mirdeep_folder}" \
        -input_linearfold_results "${linearfold_result}" \
        -output_folder "${output_folder}" 
}

produce_final_form "${FINAL_FASTA}" "${MAPPING_MIRNA_FILE}" "${DUPLICATES_INFORMATION}" "${MIRDEEP_FOLDER}" "${LINEARFOLD_RESULT}" "${OUTPUT_FOLDER}"