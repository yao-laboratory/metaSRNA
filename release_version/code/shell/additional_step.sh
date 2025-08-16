#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Global variables
CODE_PATH_AFTER_CLEAN=$1
MIN_BLOCK_LENGTH=$2
MAX_BLOCK_LENGTH=$3
MIN_BLOCK_SEQ_COUNT=$4
FORM_RESULT=$5
MAPPING_WITH_GENES_FILE=$6
OUTPUT_FOLDER=$7
echo "a$OUTPUT_FOLDER"
# Function to produce overlap analysis

# Function to check if a value is an integer
is_integer() {
    local value=$1
    [[ $value =~ ^[0-9]+$ ]]
}
if ! is_integer "$MIN_BLOCK_LENGTH" || ! is_integer "$MAX_BLOCK_LENGTH" || ! is_integer "$MIN_BLOCK_SEQ_COUNT"; then
    printf "Attention: The shortest length (%s) or longest length (%s) of the block or block sequences count threshold is not set or not a valid integer. Defaulting to 18,30,10 with order.\n" "$MIN_BLOCK_LENGTH" "$MAX_BLOCK_LENGTH" "$MIN_BLOCK_SEQ_COUNT" >&2
    MIN_BLOCK_LENGTH=18
    MAX_BLOCK_LENGTH=30
    MIN_BLOCK_SEQ_COUNT=10
    # exit 1
fi

build_additional_step() {
    local final_form=$1
    local mapping_with_genes_file=$2
    local output_folder=$3
    
    python3 "${CODE_PATH_AFTER_CLEAN}/build_additional_step.py" \
        -input_final_form "${final_form}" \
        -input_mapping_with_genes_file "${mapping_with_genes_file}" \
        -output_folder "${output_folder}" 

    python3 "${CODE_PATH_AFTER_CLEAN}/classify_mapping_patterns.py" \
        -input_bedfile "${output_folder}/sorted_representative_sequence_results.bed" \
        -min_block_length $MIN_BLOCK_LENGTH \
        -max_block_length $MAX_BLOCK_LENGTH \
        -min_block_sequences_count $MIN_BLOCK_SEQ_COUNT \
        -output_folder "${output_folder}" 
}

build_additional_step "${FORM_RESULT}" "${MAPPING_WITH_GENES_FILE}" "${OUTPUT_FOLDER}"
