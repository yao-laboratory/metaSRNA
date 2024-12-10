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

Programs and their Required Options:
  preprocess
    - Preprocess the raw input files.
    - Options:
      -r <raw_data>                       Path(Paths) to the raw data file(files, more than one file, paths use space to seperate)
      -w <preprocess functions>           Two way to choose: "merge": merge two fastq file as one fastq file; "unzip": unzip fastq.gz file to fastq file
      -o <output_folder>                  Output folder including (fastq file as extract input)
      
  extract
    - Extract data from raw input files.
    - Options:
      -r <raw_data>                       Path to the raw data file.
      -l <filter_minimum_length>          raw data fitering minimum length condition
      -o <output_folder>                  Output folder including (clean fasta file and corresponding umi file)
      -F <clean_format>                   Flexible cleaning code pattern
      --t1 <fault_torlerance>             This option sets the number of bits (errors) allowed in the pattern that are not correct.
      --t2 <tail_incomplete_tolerance>    This parameter specifies how many incomplete bits are allowed in the pattern tail.
      --umi <umi exist flag>              This flag indicates whether a UMI exists,0: not exist, 1-n:nth part in the above -F pattern

  detect_species
    - Detect Top nth species in clean fasta data.
    - Options:
      -c <clean_fasta_file> Path to the clean FASTA file.
      -t <retain how many top mapping species> Retain the top <n> mapping species' sacc and gcf values.
      -o <output_folder>    Output folder including (combined fna and (combined?) gtf files)

  detect_species_additional_step
    - download the top N species' FNA and GTF files, then combine them into combined.fna and combined.gtf.
    - Options:
      -c <gcf_files>.     File generated from detect_species step called mapping.csv.
      --cn <gcf_numbers>  Specify all GCF numbers to download, separated by commas.
      --of <output_folder>    Directory where combined.fna will be stored.
      --og <output_folder>    Directory where combined.gtf will be stored.

  map_genome
    - Map genome sequences against a reference database.
    - Options:
      -f <fna_file>          Path to the FNA file.
      -c <clean_fasta_file>  Path to the cleaned FASTA file.
      -d <database>          Path to the database directory.
      -n <database_name>     Name of the species combine database.
      -o <output_folder>     Output folder including (filter score results about mapping top species database and analysis csv about percentage of mapping)

  quantify_genome
    - Quantify the analysis results based on genome mapped data.
    - Options:
      -g <gtf_file>            Path to the GTF file.
      -m <mapping_filter_file> Path to the mapping filter score file.
      -u <umi_file>            Path to the UMI file.
      -o <output_folder>       Output folder including (quantified csv files)

  map_mirna
    - Map microRNA data to a reference.
    - Options:
      -c <clean_fasta_file> Path to the clean FASTA file.
      -d <database>         Path to the database directory.
      -n <database_name>    Name of the miRNA database.
      -o <output_folder>    Output folder including (filter score results about mapping hairpin database and analysis csv about percentage of mapping)

  quantify_mirna
    - Quantify the analysis results based on miRNA mapped data.
    - Options:
      -m <mapping_filter_file> Path to the mapping filter score file.
      -u <umi_file>            Path to the UMI file.
      -o <output_folder>       Output folder including (quantified csv files)

  integrate
    - Integrate species and miRNA data sets mapping results.
    - Options:
      -c <clean_fasta_file>  Path to the clean FASTA file.
      -s <species_file>      Path to the top species mapping score file.
      -m <mirna_file>        Path to the miRNA mapping score file.
      -u <umi_file>          Path to the UMI file.
      -o <output_folder>     Output folder including (analysis csv about how many seqid, percentage for overlapping)

  predict
    - Run predictive models on genomic data.
    - Options:
      -w <models>                choose model you wanna predict: mirdeep2, linearfold or both
      -r <raw_data>              Path to the raw genomic data.
      -c <clean_fasta_file>      Path to the clean FASTA file.
      -m <mapping_filter_score_file>  Path to the mapping filter score file.
      -f <fna_file>              Path to the FNA file.
      -n <database_name>         Name of the used database.
      -o <output_folder>         Output folder including (predictive results) 


  produce_final_form
    - Integrate species and miRNA data sets mapping results.
    - Options:
      -c  <clean_fasta_file>      Path to the clean FASTA file (after filtering length and removing duplicates)
      -m  <mapping_mirna_file>    Path to the miRNA mapping score file.
      --inf <duplicate_sequences_information> Path to file about duplicate sequences information.
      --mr <mirdeep_result>        Path to folder about mirdeep2 prediction result.
      --lf <linearfold_result>     Path to file about linearfold prediction results.
      -o <output_folder>          Output folder
    
  additional_step
    - additional_function
    - Options:
      -f  <final_form>      
      -m  <blast_score_filter_add_gene>
      -o <output_folder>          Output folder


  all
    - Execute all available processes sequentially.

Generic Options:
  -h  Display this help message

Examples:
    $0 -p extract -r ../../input/test_data1.fastq -o ../../output/extract
    $0 -p detect_species -c ../../output/extract/final_seq12.fasta -o ../../output/detect_species
    $0 -p map_genome -f ../../output/detect_species/combined_gcf.fna -c ../../output/extract/final_seq12.fa -d ../../output/map_genome/top_species_database -n species  -o ../../output/map_genome
    $0 -p map_mirna -f ../../database/microRNA_fa_database/hairpin.fa -c ../../output/extract/final_seq12.fa -d ../../output/map_mirna/hairpin_database -n hairpin  -o ../../output/map_mirna
    $0 -p integrate -s ../../output/map_genome -m ../../output/map_mirna -o ../../output/integrate
    $0 -p quantify -g ../../output/detect_species -m ../../output/map_genome/blast_score_filter.txt -u ../../output/extract/final_umi.fastq -o ../../output/quantify
    $0 -p predict -w mirdeep2 -r ../../input/test_data1.fastq  -c ../../output/extract/final_seq12.fasta -m ../../output/map_genome/blast_score_filter.txt -f ../../output/detect_species/combined_gcf.fna -n mirdeep2 -o ../../output/predict_mirdeep
    $0 -p predict -w linearfold  -c ../../output/extract/final_seq12.fasta -m ../../output/map_genome/blast_score_filter.txt  -f ../../output/detect_species/combined_gcf.fna  -o ../../output/predict_linearfold
    $0 -p predict -w both -r ../../input/test_data1.fastq  -c ../../output/extract/final_seq12.fasta -m ../../output/map_genome/blast_score_filter.txt -f ../../output/detect_species/combined_gcf.fna -n both -o ../../output/predict_mirdeep
EOF
    exit 1
}
#   $0 -p extract -f /path/to/rawdata.txt -o /path/to/output
#   $0 -p map_genome -f /path/to/genome.fna -c /path/cleaned.fasta -n dbname -o /path/to/output


# Verify and set parameters
check_params() {
    local required=("$@")
    for param in "${required[@]}"; do
        if [ -z "${opts[$param]}" ]; then
            log_error "Missing parameter: $param's value is required for $program."
            usage
        fi
    done
}

# Create directory safely
check_create_dir() {
    local dir="$1"
    local current_path=""
    IFS='/' read -r -a path_parts <<<"$dir"
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

cleanup_folder() {
    local folder="$1"

    if [ -d "$folder" ]; then
        find "$folder" -mindepth 1 -exec rm -rf {} +
        printf "Cleaned up everything in: %s\n" "$folder"
    else
        printf "The folder '%s' does not exist.\n" "$folder" >&2
    fi
}


check_create_log_dir() {
    local dir="$1"
    local program="$2"
    check_create_dir "${dir}"
    cleanup_folder "${dir}"
    check_create_dir "${dir}/log_folder"
    # Define log files with timestamps
    local timestamp=$(date +'%Y-%m-%d_%H-%M-%S')
    OUT_LOG="${dir}/log_folder/${program}_${timestamp}.out"
    ERROR_LOG="${dir}/log_folder/${program}_${timestamp}.err"
    # Initialize log files
    touch "$OUT_LOG"
    touch "$ERROR_LOG"
    log_message "Log files initialized at ${timestamp}."
}

run_script() {
    local script_path="$1"
    echo "$script_path"
    shift
    local script_name
    script_name=$(basename "$script_path")

    log_message "Starting ${script_name} step."
    if [[ -x "$script_path" ]]; then
        "$script_path" "$@" >>"$OUT_LOG" 2>>"$ERROR_LOG"
        if [[ $? -ne 0 ]]; then
            echo "Error: ${script_name} encountered an error." >&2
            return 1
        fi
    else
        echo "Error: ${script_path} does not exist or is not executable." >&2
        return 1
    fi
    log_message "Completed ${script_name} step."
}

abs_path() {
    local path
    local abs_paths
    for path in "$@"; do
        if [[ -z "$path" || "${path,,}" == "none" ]]; then
            abs_paths+="$path "
        elif [[ -e "$path" ]]; then
            local abs
            abs=$(readlink -f "$path")
            if [[ $? -ne 0 ]]; then
                echo "Error: Failed to resolve absolute path for $path" >&2
                return 1
            fi
            abs_paths+="$abs " 
        else
            echo "Error: Path does not exist - $path" >&2
            return 1
        fi
    done
    echo "${abs_paths% }"
    # echo "${abs_paths}"
}

main() {
    # Initialize variables
    local program
    declare -A opts

    # Parse command line options (with)
    while getopts ":p:f:l:o:c:n:g:u:d:m:r:s:t:w:F:h-:" opt; do
        case "$opt" in
            p) program=$OPTARG ;;
            f) opts[f]=$OPTARG ;;
            l) opts[l]=$OPTARG ;;
            o) opts[o]=$OPTARG ;;
            c) opts[c]=$OPTARG ;;
            n) opts[n]=$OPTARG ;;
            g) opts[g]=$OPTARG ;;
            u) opts[u]=$OPTARG ;;
            d) opts[d]=$OPTARG ;;
            m) opts[m]=$OPTARG ;;
            r) opts[r]=$OPTARG ;;
            s) opts[s]=$OPTARG ;;
            t) opts[t]=$OPTARG ;;
            w) opts[w]=$OPTARG ;;
            F) opts[F]=$OPTARG ;;
            -)
                case "${OPTARG}" in
                    t1)  opts[t1]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    t2)  opts[t2]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    umi) opts[umi]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    inf) opts[inf]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    mr)  opts[mr]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    lf)  opts[lf]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    cn)  opts[cn]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    of)  opts[of]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    og)  opts[og]="${!OPTIND}"; OPTIND=$((OPTIND + 1)) ;;
                    *)   log_error "Invalid option: --$OPTARG"; usage; exit 1 ;;
                esac ;;
            h) usage ;;
            \?) log_error "Invalid option: -$OPTARG"; usage ;;
            :) log_error "Option -$OPTARG requires an argument."; usage ;;
        esac
    done


    # Switch case to handle each function
    check_create_log_dir "${opts[o]}" "$program"
    
    local cleaning_code_path=${PARENT_DIR}/python/cleaning
    local after_cleaning_code_path=${PARENT_DIR}/python/after_cleaning
    local refprok_database=${RELEASE_DIR}/database/blast_refprok_database
    local param_flag=0

    case "$program" in
        preprocess)
            check_params r w o
            check_create_dir "${opts[o]}/middle_results"
            IFS=',' read -r -a input_paths <<< "${opts[r]}"
            printf "path_r: %s\n" "${opts[r]}"
            printf "path: %s\n" "${input_paths[@]}"
            expanded_paths=()
            for path in "${input_paths[@]}"; do
                abs_path_result=$(abs_path "$path")
                printf "abs_path_result: %s\n" "$abs_path_result"
                expanded_paths+=("$abs_path_result")
            done
            run_script "${DIR}/preprocess.sh" "$cleaning_code_path" "${opts[w]}" "$(abs_path "${opts[o]}")" "${#expanded_paths[@]}" "${expanded_paths[@]}"
            ;;

        extract)
            check_params r l o F
            if [ "${opts[F]}" != "clean" ]; then
                check_params t1 t2 umi
                t1_param="${opts[t1]}"
                t2_param="${opts[t2]}"
                umi_param="${opts[umi]}"
            else
                t1_param=""
                t2_param=""
                umi_param="0"
            fi
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[r]}" "${opts[o]}")"
            run_script "${DIR}/run.sh" "$cleaning_code_path" "${opts[l]}" "${opts[F]}" "$t1_param" "$t2_param" "$umi_param" "${paths[@]}"
            ;;

        detect_species)
            check_params c o t
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[c]}" "${opts[o]}")"
            run_script "${DIR}/species_detection.sh" "$refprok_database" "$after_cleaning_code_path" "${opts[t]}" "${paths[@]}"
            ;;
        
        detect_species_additional_step)
            check_params o of og
            check_create_dir "${opts[o]}/middle_results"
            if [[ -n "${opts[c]}" ]]; then
                param_flag=0
                IFS=' ' read -r -a paths <<< "$(abs_path "${opts[c]}" "${opts[o]}" "${opts[of]}" "${opts[og]}")"
                run_script "${DIR}/species_detection_additional_step.sh" "$param_flag" "${paths[@]}"
            else
                check_params cn
                param_flag=1
                IFS=' ' read -r -a paths <<< "$(abs_path "${opts[o]}" "${opts[of]}" "${opts[og]}")"
                run_script "${DIR}/species_detection_additional_step.sh" "$param_flag" "${opts[cn]}" "${paths[@]}"
            fi
            ;;
        
        map_genome)
            check_params f c d n o
            check_create_dir "${opts[d]}"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[f]}" "${opts[c]}" "${opts[d]}" "${opts[o]}")"
            run_script "${DIR}/map_genome.sh" "${opts[n]}" "${paths[@]}"
            ;;

        quantify_genome)
            check_params g m u o
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[g]}" "${opts[m]}" "${opts[u]}" "${opts[o]}")"
            run_script "${DIR}/species_function_quantification.sh" "$after_cleaning_code_path" "${paths[@]}"
            
            ;;

        map_mirna)
            check_params f c d n o
            check_create_dir "${opts[d]}"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[f]}" "${opts[c]}" "${opts[d]}" "${opts[o]}")"
            run_script "${DIR}/map_microrna.sh" "${opts[n]}" "$after_cleaning_code_path" "${paths[@]}" 

            ;;

        quantify_mirna)
            check_params m u o
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[m]}" "${opts[u]}" "${opts[o]}")"
            run_script "${DIR}/mirna_function_quantification.sh" "$after_cleaning_code_path" "${paths[@]}"

            ;;
        integrate)
            check_params c s m u o
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[c]}" "${opts[s]}" "${opts[m]}" "${opts[u]}" "${opts[o]}")"
            run_script "${DIR}/integration.sh" "$cleaning_code_path" "$after_cleaning_code_path" "${paths[@]}" 
            ;;

        predict)
            check_params w
            check_create_dir "${opts[o]}/middle_results"
            check_create_dir "${opts[o]}/middle_results/database"
            
            run_mirdeep2() {
                check_params c f n o
                IFS=' ' read -r -a paths <<< "$(abs_path "${opts[c]}" "${opts[f]}" "${opts[o]}")"
                run_script "${DIR}/run_mirdeep2.sh" "${opts[n]}" "$after_cleaning_code_path" "${paths[@]}" 
            }

            run_linearfold() {
                check_params c m f o
                IFS=' ' read -r -a paths <<< "$(abs_path "${opts[c]}" "${opts[m]}" "${opts[f]}" "${opts[o]}")"
                run_script "${DIR}/run_linearfold.sh"  "$after_cleaning_code_path" "${paths[@]}" 
            }

            model=${opts[w]}
            case "$model" in
                mirdeep2)
                    run_mirdeep2
                    ;;
                linearfold)
                    run_linearfold
                    ;;
                both)
                    log_message "Starting predict steps sequentially."
                    run_mirdeep2
                    run_linearfold
                    log_message "Completed all predict steps."
                    ;;
                *)
                    printf "Error: Invalid value for choosing models parameter: %s\n" "$model" >&2
                    return 1
                    ;;
            esac
            ;;

        produce_final_form)
            check_params c m inf mr lf o
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[c]}" "${opts[m]}" "${opts[inf]}" "${opts[mr]}" "${opts[lf]}" "${opts[o]}")"
            run_script "${DIR}/final_form_production.sh" "$after_cleaning_code_path" "${paths[@]}"
            ;;
        
        additional_step)
            check_params f m
            check_create_dir "${opts[o]}/middle_results"
            IFS=' ' read -r -a paths <<< "$(abs_path "${opts[f]}" "${opts[m]}" "${opts[o]}")"
            run_script "${DIR}/additional_step.sh" "$after_cleaning_code_path" "${paths[@]}"
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
}

main "$@"