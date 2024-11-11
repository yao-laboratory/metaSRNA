#!/bin/bash
# Exit immediately if any command fails
set -e

CODE_PATH=$1
FUNCTION=$2
OUTPUT_FOLDER=$3
FILES_NUM=$4
shift 4  # Remove the first four known parameters

# Store the paths into an array for further use
paths=()

if [[ "$FUNCTION" == "merge" ]]; then
    if [[ "$FILES_NUM" -eq 2 ]]; then
        file_1=$1
        file_2=$2
        shift 2

        # # strip the path and keep only the file name
        # fastq_1=$(basename "$file_1")
        # fastq_2=$(basename "$file_2")

        # if [[ "$fastq_1" == *"_"* ]]; then
        #     merge_output="${OUTPUT_FOLDER}/${fastq_1%_*}_pair.fastq"
        # else
        #     merge_output="${OUTPUT_FOLDER}/${fastq_1}_pair.fastq"
        # fi

        bbmerge.sh in1="$file_1" in2="$file_2" out="${OUTPUT_FOLDER}/merged.fastq"
    else
        printf "Error: The 'merge' function requires exactly 2 input files.\n" >&2
        exit 1
    fi

elif [[ "$FUNCTION" == "unzip" ]]; then
    for ((i = 0; i < FILES_NUM; i++)); do
        paths+=("$1")
        shift 1
    done

    printf "Number of paths provided: %d\n" "$FILES_NUM"

    for path in "${paths[@]}"; do
        file=$(basename "$path")
        output_path="${OUTPUT_FOLDER}/${file%.fastq.gz}.fastq"
        printf "Unzipping: %s -> %s\n" "$path" "$output_path"
        gzip -d -c "$path" > "$output_path"
    done
else
    printf "Error: Unsupported function '%s'.\n" "$FUNCTION" >&2
    exit 1
fi



