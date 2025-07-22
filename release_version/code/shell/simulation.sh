#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

# Global variables
CODE_PATH_AFTER_CLEAN=$1
BLOCKS_NUM=$2
MIN_GAP=$3
MAX_GAP=$4
SIMULATION_TIMES=$5
DRAW_FLAG=$6
OUTPUT_FOLDER=$7
FILES_NUM=$8
shift 8  
INPUT_FILES="$*"
echo "a$INPUT_FILES"
echo "a$OUTPUT_FOLDER"
# Function to produce overlap analysis
simulate_blocks() {
    local input_files=$1
    local blocks_number=$2
    local min_blocks_gap=$3
    local max_blocks_gap=$4
    local simulation_times=$5
    local draw_flag=$6
    local output_folder=$7
    
    python3 "${CODE_PATH_AFTER_CLEAN}/simulation.py" \
        -input_files "${input_files}" \
        -number_of_blocks "${blocks_number}" \
        -block_gap_min_threshold "${min_blocks_gap}" \
        -block_gap_max_threshold "${max_blocks_gap}" \
        -simulation_times "${simulation_times}" \
        -drawing_flag "${draw_flag}" \
        -output_folder "${output_folder}" \
        -code_folder "${CODE_PATH_AFTER_CLEAN}"

}

simulate_blocks "${INPUT_FILES}" "${BLOCKS_NUM}" "${MIN_GAP}" "${MAX_GAP}" "${SIMULATION_TIMES}" "${DRAW_FLAG}" "${OUTPUT_FOLDER}"
