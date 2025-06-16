#!/bin/bash

# Exit immediately if any command fails
set -e
set -o pipefail
# module load blast/2.14
# module load biodata/1.0
# Global variables for input arguments
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
    printf "Error: FLAG must be a numeric value.\n" >&2
    exit 1
fi
FLAG="$1"
RESULTS="$3"
FNA_FOLDER="$4"
GTF_FOLDER="$5"
FILE_NAME="$6"
MIDDLE_FOLDER="$RESULTS/middle_results"
# mkdir -p "$MIDDLE_FOLDER"mkdir -p "$MIDDLE_FOLDER"
BLAST_REFPROK_MAPPING="$MIDDLE_FOLDER/ncbi_dataset/ncbi_dataset/data"

get_gcfnumbers() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        printf "Error: File not found: %s\n" "$file" >&2
        return 1
    fi

    local gcf_string
    gcf_string=$(awk -F'\t' 'NR > 1 { printf "%s ", $2 }' "$file" | sed 's/ $//')

    if [[ -z "$gcf_string" ]]; then
        printf "Error: No GCF values found in the file.\n" >&2
        return 1
    fi

    printf "%s\n" "$gcf_string"
}
# Function to download genomes using accession numbers from mapping.csv
download_genomes() {
  local gcf_numbers="$1"
  gcf_numbers=$(printf "%s" "$GCF_NUMBERS" | tr ',' ' ')
  # gcf_numbers=$(awk -F'\t' 'NR > 1 && $2 {print $2}' "$RESULTS/mapping.csv" | tr '\n' ' ')
  
  if [[ -z "$gcf_numbers" ]]; then
    printf "No GCF numbers found. Exiting.\n" >&2
    return 1
  fi
  
  printf "GCF numbers: %s\n" "$gcf_numbers"

  # Ensure dataset CLI is downloaded
  if [[ ! -f "$MIDDLE_FOLDER/datasets" ]]; then
    curl -o "$MIDDLE_FOLDER/datasets" 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
    chmod 777 "$MIDDLE_FOLDER/datasets"
  fi

  "${MIDDLE_FOLDER}/datasets" download genome accession ${gcf_numbers} \
    --include genome,gtf --filename "$MIDDLE_FOLDER/ncbi_dataset.zip"

  unzip -o "$MIDDLE_FOLDER/ncbi_dataset.zip" -d "$MIDDLE_FOLDER/ncbi_dataset"

  find "$BLAST_REFPROK_MAPPING" -type f -name "*.fna" -exec cat {} + > "$RESULTS/combined_${FILE_NAME}.fna"
  find "$BLAST_REFPROK_MAPPING" -type f -name "*.gtf" -exec cat {} + > "$RESULTS/combined_${FILE_NAME}.gtf"


  # Copy files to respective folders
  printf "Copying fna to '%s'...\n" "$FNA_FOLDER"
  cp -f "$RESULTS/combined_${FILE_NAME}.fna" "$FNA_FOLDER" || {
      printf "Error: Failed to copy fna file to '%s'.\n" "$FNA_FOLDER" >&2
      exit 1
  }

  printf "Copying gtf to '%s'...\n" "$GTF_FOLDER"
  cp -f "$RESULTS/combined_${FILE_NAME}.gtf" "$GTF_FOLDER" || {
      printf "Error: Failed to copy gtf file to '%s'.\n" "$GTF_FOLDER" >&2
      exit 1
  }
  printf "species detect addtional step finished.\n"
}

# Main function to coordinate the script flow
main() {
    local GCF_NUMBERS GCF_FILE

    if (( $FLAG == 1 )); then
        GCF_NUMBERS="$2"
    elif (( $FLAG == 0 )); then
        GCF_FILE="$2"
        GCF_NUMBERS=$(get_gcfnumbers "$GCF_FILE")
        if [[ -z "$GCF_NUMBERS" ]]; then
            printf "Error: Failed to get GCF numbers from file %s\n" "$GCF_FILE" >&2
            return 1
        fi
    else
        printf "Error: Invalid FLAG value: %d\n" "$FLAG" >&2
        return 1
    fi

    download_genomes "$GCF_NUMBERS"
}

# Execute main function
main "$@"
