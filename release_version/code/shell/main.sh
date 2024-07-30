#!/bin/bash

# # Activate the environment
# source activate cleanfasq_env

# Define this shell directory and its' parent directory
DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
PARENT_DIR="$(dirname "$DIR")"
RELEASE_DIR="$(dirname "$PARENT_DIR")"

echo "Running main assembly process in the folder $DIR"

# Define two log fileS
OUT_LOG=""
ERROR_LOG=""

# Function to log messages and errors
log_message() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1" | tee -a "$OUT_LOG"
}

log_error() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ERROR: $1" | tee -a "$ERROR_LOG"
}

# Function to display help/usage information
usage() {
    cat <<EOF
Usage: $0 -p <program> [options]
EOF
    exit 1
}

# Initialize variables
program=""
declare -A opts

# Parse command line options
while getopts ":p:f:o:c:n:g:u:d:m:r:s:h" opt; do
    case "$opt" in
        p) program=$OPTARG ;;
        f) opts[f]=$OPTARG ;;
        o) opts[o]=$OPTARG ;;
        c) opts[c]=$OPTARG ;;
        n) opts[n]=$OPTARG ;;
        g) opts[g]=$OPTARG ;;
        u) opts[u]=$OPTARG ;;
        d) opts[d]=$OPTARG ;;
        m) opts[m]=$OPTARG ;;
        r) opts[r]=$OPTARG ;;
        s) opts[s]=$OPTARG ;;
        h) usage ;;
    esac
done

# Verify and set parameters based on program
check_params() {
    local required=("$@")
    for param in "${required[@]}"; do
        if [ -z "${opts[$param]}" ]; then
            log_error "Missing parameter: $param is required for $program."
            usage
        fi
    done
}

# Create directory safely
check_create_dir() {
    local dir="$1"
    local current_path=""
    IFS='/' read -r -a path_parts <<< "$dir"
    for part in "${path_parts[@]}"; do
        if [[ -z "$current_path" ]]; then
            current_path="$part"
        else
            current_path="$current_path/$part"
        fi

        if [[ ! -d "$current_path" ]]; then
            if mkdir -p "$current_path"; then
                printf "Directory created: %s\n" "$current_path"
            fi
        fi
    done
}

check_create_log_dir() {
    local dir="$1"
    local program="$2"
    check_create_dir "$dir"
    # Define log files with timestamps
    local timestamp=$(date +'%Y-%m-%d_%H-%M-%S')
    OUT_LOG="${dir}/${program}_${timestamp}.out"
    ERROR_LOG="${dir}/${program}_${timestamp}.err"

    # Initialize log files
    touch "$OUT_LOG"
    touch "$ERROR_LOG"
    log_message "Log files initialized at ${timestamp}."
}


# Switch case to handle each program
case "$program" in
    extract)
        check_params r o
        check_create_log_dir "${opts[o]}/log_folder" "$program"
        check_create_dir "${opts[o]}/middle_results"
        log_message "Starting extract step."
        if [[ -x "${DIR}/run.sh" ]]; then
            "${DIR}/run.sh" "$(readlink -f "${opts[r]}")" "$(readlink -f "${opts[o]}")" "${PARENT_DIR}/python/cleaning" >> "$OUT_LOG" 2>> "$ERROR_LOG"
        else
            echo "Error: ${DIR}/run.sh does not exist or is not executable." >&2
            exit 1
        fi
        log_message "Completed extract step."
        ;;

    detect_species)
        check_params c o
        check_create_log_dir "${opts[o]}/log_folder" "$program"
        check_create_dir "${opts[o]}/middle_results"
        log_message "Starting extract step."
        if [[ -x "${DIR}/species_detection.sh" ]]; then
            echo "${DATABASE_DIR}/database/blast_refprok_database"
            "${DIR}/species_detection.sh" "$(readlink -f "${opts[c]}")" "$(readlink -f "${opts[o]}")" "${RELEASE_DIR}/database/blast_refprok_database" "${PARENT_DIR}/python/after_cleaning">> "$OUT_LOG" 2>> "$ERROR_LOG"
        else
            echo "Error: ${DIR}/species_detection.sh does not exist or is not executable." >&2
            exit 1
        fi
        log_message "Completed detect_species step."
        ;;
        
    map_genome)
        check_params f c d n o
        check_create_log_dir "${opts[o]}/log_folder"
        check_create_dir "${opts[d]}"
        log_message "Starting map_genome step."
        if [[ -x "${DIR}/genome_mapping.sh" ]]; then
            # echo "${DIR}/genome_mapping.sh"
            # echo "$(readlink -f "${opts[f]}")"  
            # echo "$(readlink -f "${opts[c]}")"  
            # echo "$(readlink -f "${opts[d]}")"  
            # echo "${opts[n]}" 
            # echo "$(readlink -f "${opts[o]}")"
            "${DIR}/genome_mapping.sh" "$(readlink -f "${opts[f]}")"  "$(readlink -f "${opts[c]}")"  "$(readlink -f "${opts[d]}")"  "${opts[n]}" "$(readlink -f "${opts[o]}")" >> "$OUT_LOG" 2>> "$ERROR_LOG"
        else
            echo "Error: ${DIR}/genome_mapping.sh does not exist or is not executable." >&2
            exit 1
        fi
        log_message "Completed map_genome step."
        ;;
        
    quantify)
        check_params g m u o
        check_create_log_dir "${opts[o]}/log_folder"
        check_create_dir "${opts[o]}/middle_results"
        log_message "Starting quantify step."
        if [[ -x "${DIR}/function_quantification.sh" ]]; then
            "${DIR}/function_quantification.sh" "$(readlink -f "${opts[g]}")" "$(readlink -f "${opts[m]}")" "$(readlink -f "${opts[u]}")" "$(readlink -f "${opts[o]}")" "${PARENT_DIR}/python/after_cleaning" >> "$OUT_LOG" 2>> "$ERROR_LOG"
        else
            echo "Error: ${DIR}/function_quantification.sh does not exist or is not executable." >&2
            exit 1
        fi
        log_message "Completed quantify step."
        ;;
        
    map_mirna)
        check_params f c d n o
        check_create_log_dir "${opts[o]}/log_folder"
        check_create_dir "${opts[d]}"
        log_message "Starting map_miRNA step."
        if [[ -x "${DIR}/map_microrna.sh" ]]; then
            "${DIR}/map_microrna.sh" "$(readlink -f "${opts[f]}")"  "$(readlink -f "${opts[c]}")"  "$(readlink -f "${opts[d]}")"  "${opts[n]}"  "${PARENT_DIR}/python/after_cleaning" "$(readlink -f "${opts[o]}")" >> "$OUT_LOG" 2>> "$ERROR_LOG"
        else
            echo "Error: ${DIR}/map_microrna.sh does not exist or is not executable." >&2
            exit 1
        fi
        log_message "Completed map_miRNA step."
        ;;
        
    integrate)
        check_params s m o
        check_create_log_dir "${opts[o]}/log_folder"
        log_message "Starting integrate step."
        if [[ -x "${DIR}/integration.sh" ]]; then
            "${DIR}/integration.sh" "$(readlink -f "${opts[s]}")" "$(readlink -f "${opts[m]}")" "${PARENT_DIR}/python/after_cleaning" "$(readlink -f "${opts[o]}")" >> "$OUT_LOG" 2>> "$ERROR_LOG"
        else
            echo "Error: ${DIR}/integration.sh does not exist or is not executable." >&2
            exit 1
        fi
        log_message "Completed integrate step."
        ;;
        
    predict)
        check_params r c m f n o
        check_create_log_dir "${opts[o]}/log_folder"
        check_create_dir "${opts[o]}/middle_results"
        log_message "Starting predict steps simultaneously."
        "${DIR}/run_mirdeep2.sh" "${opts[c]}" "${opts[m]}" "${opts[f]}" "${opts[n]}" "${PARENT_DIR}/code/python/after_cleaning" "${opts[o]}/middle_results" "${opts[o]}" >> "$OUT_LOG" 2>> "$ERR_LOG" &
        mirdeep2_pid=$!
        "${DIR}/run_linearfold.sh" "${opts[c]}" "${opts[m]}" "${PARENT_DIR}/code/python/after_cleaning" "${opts[o]}/middle_results" "${opts[o]}" >> "$OUT_LOG" 2>> "$ERR_LOG" &
        linearfold_pid=$!
        wait $mirdeep2_pid
        wait $linearfold_pid
        log_message "Completed predict steps."
        ;;
        
    all)
        log_message "Running all processes."
        # Implement calls to all functions here
        ;;
        
    *)
        log_error "Unknown command: '$program'."
        usage
        ;;
esac
