#!/bin/bash
# Exit immediately if any command fails
set -euo pipefail

# Global variables
CODE_PATH=$1
GTF_FILE=$2
FILTER_SCORE_FILE=$3
UMI_FILE=$4
RESULTS_DIR=$5
MIDDLE_RESULTS_DIR=${RESULTS_DIR}/middle_results

filter_gtf_comments() {
    local input_file=$1
    local temp_file

    # Create a secure temporary file
    temp_file=$(mktemp)

    # Use grep to filter out lines starting with '#' and write to the temporary file
    if ! grep -v '^#' "$input_file" > "$temp_file"; then
        echo "Error: Failed to filter comments from $input_file" >&2
        rm -f "$temp_file"  # Clean up temporary file in case of failure
        return 1
    fi

    # Move the temporary file back to the original file
    if ! mv "$temp_file" "$input_file"; then
        echo "Error: Failed to overwrite the original file $input_file" >&2
        rm -f "$temp_file"  # Clean up temporary file in case of failure
        return 1
    fi
}

run_blast_process() {
    local operation=$1
    shift
    local args=("$@")
    
    if ! python3 "${CODE_PATH}/blast_process_main.py" "$operation" "${args[@]}"; then
        printf "Error: Failed to run blast_process_main.py %s\n" "$operation" >&2
        return 1
    fi
}

quantify_species(){
    local sorted_score_file=$1
    declare -A unique_output_files=(
            ["gene_biotype"]="${RESULTS_DIR}/output_unique_biotype.csv"
            ["gene gene_biotype"]="${RESULTS_DIR}/output_unique_gene_biotype.csv"
            ["gene"]="${RESULTS_DIR}/output_unique_gene.csv"
            )

    for keys in "${!unique_output_files[@]}"; do
        run_blast_process "find_unique_count" \
            -input_sorted_score_file "$sorted_score_file" \
            -input_umi_file "$UMI_FILE" \
            -input_keys "$keys" \
            -output_unique "${unique_output_files[$keys]}"
    done
}

# Main function
main() {

    # filter GTF file to remove comments
    filter_gtf_comments "$GTF_FILE"

    # add genes information to blast results
    run_blast_process "add_gene" \
        -input_gcf "$GTF_FILE" \
        -input_score "$FILTER_SCORE_FILE" \
        -output_csv "${MIDDLE_RESULTS_DIR}/blast_score_filter_add_gene.csv"
    # quantify species
    quantify_species "${MIDDLE_RESULTS_DIR}/blast_score_filter_add_gene.csv"

   
}

# Execute main function
main "$@"