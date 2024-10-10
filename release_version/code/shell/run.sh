#!/bin/bash
# Exit immediately if any command fails
set -e

CODE_PATH=$1
NUM=$2
FORMAT=$3
FAULT_TL=$4
INCOMPLETE_TL=$5
UMI_FLAG=$6
INPUT=$7
RESULTS=$8

# echo "$FORMAT"
# echo "$FAULT_TL"
# echo "$INCOMPLETE_TL"
echo "$INPUT"
# echo "$RESULTS"
# echo "$UMI_FLAG"
# echo "$CODE_PATH"

# run clean_command script
run_clean_command() {
    if [[ "$FORMAT" != "clean" ]]; then
        python3 "${CODE_PATH}/clean_command.py" clean_fastq \
            -input "$INPUT" \
            -output_filename "${RESULTS}/middle_results/output" \
            -fa_format "$FORMAT" \
            -fault_tolerance "$FAULT_TL" \
            -tail_incomplete_tolerance "$INCOMPLETE_TL"
    else
        python3 "${CODE_PATH}/clean_command.py" after_clean \
            -input "$INPUT" \
            -output_path "${RESULTS}/middle_results/final_seq.fastq"
    fi
}

# produce original sequence file
produce_seq() {
    local step1="${RESULTS}/middle_results/output_1_step1.fastq"
    local step2="${RESULTS}/middle_results/output_1_step2.fastq"
    local final_seq="${RESULTS}/middle_results/final_seq.fastq"
    local output_folder="$RESULTS"
    if [[ "$FORMAT" != "clean" ]]; then
        if [[ -f "$step2" ]]; then
            cat "$step1" "$step2" > "$final_seq"
            printf "Concatenated '%s' and '%s' into '%s'.\n" "$step1" "$step2" "$final_seq"
        else
            cat "$step1" > "$final_seq"
            printf "File '%s' does not exist, concatenating only '%s' into '%s' due to fault tolerance is 0.\n" "$step2" "$step1" "$final_seq"
        fi
    fi

    python3 "${CODE_PATH}/filter_sequence_length.py" process_length \
        -input "$final_seq" \
        -output_folder "$output_folder" \
        -filter_min_length "$NUM"

    # Convert FASTQ to FASTA formats
    seqtk seq -a "${RESULTS}/final_seq_${NUM}.fastq" > "${RESULTS}/final_seq_${NUM}.fasta"
    seqtk seq -a "${RESULTS}/final_seq_${NUM}.fastq" > "${RESULTS}/final_seq_${NUM}.fa"
}


# produce umi file
produce_umi() {
    if [[ "$UMI_FLAG" != 0 ]]; then
        local umi_step1="${RESULTS}/middle_results/output_${UMI_FLAG}_step1.fastq"
        local umi_step2="${RESULTS}/middle_results/output_${UMI_FLAG}_step2.fastq"
        local final_umi="${RESULTS}/final_umi.fastq"
        local final_umi_fasta="${RESULTS}/final_umi.fasta"

        cat "$umi_step1" "$umi_step2" > "$final_umi"
        seqtk seq -a "$final_umi" > "$final_umi_fasta"
        
        # De-duplication: same seq and same umi
        # python3 "${CODE_PATH}/clean_command.py" De-duplication ..
    fi
}
main() {
    run_clean_command
    produce_seq
    produce_umi
    # wc -l $results/$1.fastq | awk '{print $1/4}'
    # wc -l $results/final.fastq | awk '{print $1/4}'
    # wc -l $results/final_seq12.fastq | awk '{print $1/4}'
}
main "$@"

