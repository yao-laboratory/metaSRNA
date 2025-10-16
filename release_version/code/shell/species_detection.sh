#!/bin/bash

# Exit immediately if any command fails
set -e
set -o pipefail

# Global variables for input arguments

CODE_PATH="$1"
TOP_NUM="$2"
CLEAN_FASTA="$3"
RESULTS="$4"
REFPROK_DATABASE="$5"

MIDDLE_FOLDER="$RESULTS/middle_results"

printf "Refprok Database: %s\nCode Path: %s\nClean Fasta: %s\nResults: %s\n" \
  "$REFPROK_DATABASE" "$CODE_PATH" "$CLEAN_FASTA" "$RESULTS"

preprocess_random_select() {
  local input_file="$CLEAN_FASTA"
  local output_file="${MIDDLE_FOLDER}/random_selected_final_seq_12.fasta"

  python3 "${CODE_PATH}/preprocess_random_select.py" \
    -input_file "$input_file" \
    -output_file "$output_file" \
    -select_proportion_divisor 200
}

# run blastn
run_blastn() {
  local blast_output="$MIDDLE_FOLDER/blastn_refprok_evalue10.txt"
  local seleted_clean_fasta="${MIDDLE_FOLDER}/random_selected_final_seq_12.fasta"
  blastn -query "$seleted_clean_fasta" -db "${REFPROK_DATABASE}/ref_prok_rep_genomes" \
    -out "$blast_output" -num_threads 4 -evalue 10 \
    -outfmt "6 qseqid sacc sstart send evalue bitscore qcovhsp pident sseqid staxids sscinames sblastnames" \
    -max_target_seqs 5 -max_hsps 1 -task "blastn"
  
  if [[ ! -s "$blast_output" ]]; then
    printf "BLASTn output is empty. Exiting.\n" >&2
    return 1
  fi

  printf "BLASTn complete. Output written to %s\n" "$blast_output"
}

# Function to perform line count and unique sequence extraction
process_blast_output() {
  local blast_output="$MIDDLE_FOLDER/blastn_refprok_evalue10.txt"
  local unique_seq_count

  printf "Processing BLAST output file: %s\n" "$blast_output"
  wc -l "$blast_output"
  unique_seq_count=$(awk -F'\t' '{print $1}' "$blast_output" | sort -n | uniq -c | wc -l)
  
  printf "Unique sequence count: %s\n" "$unique_seq_count"
}

# Function to run species classification scripts
run_species_classification() {
  local num_sequences
  local species_output
  local seleted_clean_fasta="${MIDDLE_FOLDER}/random_selected_final_seq_12.fasta"
  num_sequences=$(wc -l $seleted_clean_fasta |awk '{print $1/2}')
  printf "Total sequences: %s\n" "$num_sequences"

  python3 "${CODE_PATH}/blast_species_classification.py" \
    -input_folder "$MIDDLE_FOLDER" \
    -total_sequence_number "$num_sequences" \
    -output_folder "$RESULTS"

  python3 "${CODE_PATH}/union_top_species_classification.py" \
    -input_folder "$MIDDLE_FOLDER" \
    -top_count "$TOP_NUM" \
    -output_folder "$RESULTS"
}


# Main function to coordinate the script flow
main() {
  preprocess_random_select
  run_blastn
  process_blast_output
  run_species_classification
}

# Execute main function
main "$@"
