#!/bin/bash
# Exit immediately if any command fails
#!/bin/bash

set -euo pipefail

# Global variables
CODE_PATH=$1
FILTER_SCORE_FILE=$2
UMI_FILE=$3
RESULTS_DIR=$4
MIDDLE_RESULTS_DIR=${RESULTS_DIR}/middle_results

copy_score_file() {
    local source_file="$1"
    local destination_file="$2"
    
    if ! cp "$source_file" "$destination_file"; then
        printf "Error: Failed to copy %s to %s\n" "$source_file" "$destination_file" >&2
        return 1
    fi
}

run_mirna_quantification() {
    local input_score_file="$1"
    local umi_file="$2"
    local output_file="$3"
    
    if ! python3 "${CODE_PATH}/blast_process_main.py" find_unique_count_mirna \
        -input_sorted_score_file "$input_score_file" \
        -input_umi_file "$umi_file" \
        -output_unique "$output_file"; then
        printf "Error: miRNA quantification process failed.\n" >&2
        return 1
    fi
}

main() {
    local blastn_hairpinrna_file="${MIDDLE_RESULTS_DIR}/blastn_hairpinrna.txt"
    local output_unique_biotype_mirna="${RESULTS_DIR}/output_unique_biotype_mirna.csv"

    # copy file
    copy_score_file "$FILTER_SCORE_FILE" "$blastn_hairpinrna_file"

    # run miRNA quantification
    run_mirna_quantification "$blastn_hairpinrna_file" "$UMI_FILE" "$output_unique_biotype_mirna"
}

main "$@"
