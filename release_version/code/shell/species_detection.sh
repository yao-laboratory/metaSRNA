#!/bin/bash

# Exit immediately if any command fails
set -e
set -o pipefail
module load blast/2.14
module load biodata/1.0
# Global variables for input arguments
REFPROK_DATABASE="$1"
CODE_PATH="$2"
TOP_NUM="$3"
CLEAN_FASTA="$4"
RESULTS="$5"

MIDDLE_FOLDER="$RESULTS/middle_results"
# mkdir -p "$MIDDLE_FOLDER"mkdir -p "$MIDDLE_FOLDER"
BLAST_REFPROK_MAPPING_FNA="$MIDDLE_FOLDER/ncbi_dataset_fna/ncbi_dataset/data"

printf "Refprok Database: %s\nCode Path: %s\nClean Fasta: %s\nResults: %s\n" \
  "$REFPROK_DATABASE" "$CODE_PATH" "$CLEAN_FASTA" "$RESULTS"

# run blastn
run_blastn() {
  local blast_output="$MIDDLE_FOLDER/blastn_refprok_evalue10.txt"
  blastn -query "$CLEAN_FASTA" -db "${REFPROK_DATABASE}/ref_prok_rep_genomes" \
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

  num_sequences=$(wc -l $CLEAN_FASTA |awk '{print $1/2}')
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

# Function to download genomes using accession numbers from mapping.csv
download_genomes() {
  local gcf_numbers

  gcf_numbers=$(awk -F'\t' 'NR > 1 && $2 {print $2}' "$RESULTS/mapping.csv" | tr '\n' ' ')
  
  if [[ -z "$gcf_numbers" ]]; then
    printf "No GCF numbers found in mapping.csv. Exiting.\n" >&2
    return 1
  fi
  
  printf "GCF numbers: %s\n" "$gcf_numbers"

  # Ensure dataset CLI is downloaded
  if [[ ! -f "$MIDDLE_FOLDER/datasets" ]]; then
    curl -o "$MIDDLE_FOLDER/datasets" 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
    chmod 777 "$MIDDLE_FOLDER/datasets"
  fi

  "${MIDDLE_FOLDER}/datasets" download genome accession ${gcf_numbers} \
    --include genome --filename "$MIDDLE_FOLDER/ncbi_dataset_fna.zip"

  unzip -o "$MIDDLE_FOLDER/ncbi_dataset_fna.zip" -d "$MIDDLE_FOLDER/ncbi_dataset_fna"

  find "$BLAST_REFPROK_MAPPING_FNA" -type f -name "*.fna" -exec cat {} + > "$RESULTS/combined_gcf.fna"
}

# Main function to coordinate the script flow
main() {
  run_blastn
  process_blast_output
  run_species_classification
  download_genomes
}

# Execute main function
main "$@"
