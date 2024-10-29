#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Global variables
CODE_PATH_AFTER_CLEAN=$1
FORM_RESULT=$2
MAPPING_WITH_GENES_FILE=$3
OUTPUT_FOLDER=$4
echo "a$OUTPUT_FOLDER"
# Function to produce overlap analysis
build_additional_step() {
    local final_form=$1
    local mapping_with_genes_file=$2
    local output_folder=$3
    
    python3 "${CODE_PATH_AFTER_CLEAN}/build_additional_step.py" \
        -input_final_form "${final_form}" \
        -input_mapping_with_genes_file "${mapping_with_genes_file}" \
        -output_folder "${output_folder}" 
}

build_additional_step "${FORM_RESULT}" "${MAPPING_WITH_GENES_FILE}" "${OUTPUT_FOLDER}"