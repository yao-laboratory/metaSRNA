import argparse
import time
from add_gene_information import add_gene
from find_seq_umi_count import find_unique_count
from find_seq_umi_count import find_unique_count_mirna
from extract_sequences_from_genome import extract_sequences_from_genome
from filter_mapping_result import keep_best_blast_hits
    
def main():
    parser = argparse.ArgumentParser(prog='blastprocess')

    subparsers = parser.add_subparsers(
        dest='subcommand', help='Sub Commands Help')
    # add sub command
    parser_c1 = subparsers.add_parser("add_gene",
                                      help='')
    parser_c1.add_argument('-input_gcf', required=True,
                            type=str, help='input file: gtf', default="none")
    parser_c1.add_argument('-input_score', required=True,
                            type=str, help='input file: txt', default="none")
    parser_c1.add_argument('-temp_folder', required=True,
                            type=str, help='temp files folder', default="none")
    parser_c1.add_argument('-output_csv', required=True,
                            type=str, help='output file: csv', default="none")
    
    parser_c2 = subparsers.add_parser("find_unique_count",
                                      help='')
    
    parser_c2.add_argument('-input_sorted_score_file', required=True,
                            type=str, help='input file: add_gene parser output csv format', default="none")
    
    parser_c2.add_argument('-input_umi_file', required=True,
                            type=str, help='input file: fastq format', default="none")

    parser_c2.add_argument('-input_keys', required=True,
                            type=str, help='', default="none")

    parser_c2.add_argument('-output_unique', required=True,
                            type=str, help='output file: csv format', default="none")
    
    parser_c3 = subparsers.add_parser("find_unique_count_mirna",
                                      help='')
    
    parser_c3.add_argument('-input_sorted_score_file', required=True,
                            type=str, help='input file: add_gene parser output csv format', default="none")
    
    parser_c3.add_argument('-input_umi_file', required=True,
                            type=str, help='input file: fastq format', default="none")

    parser_c3.add_argument('-output_unique', required=True,
                            type=str, help='output file: csv format', default="none")

    parser_c4 = subparsers.add_parser("extract_sequences_from_genome",
                                      help='')
    
    parser_c4.add_argument('-input_blast_result', required=True,
                            type=str, help='input blast filter file', default="none")
    
    parser_c4.add_argument('-input_genome', required=True,
                            type=str, help='input fasta file', default="none")

    parser_c4.add_argument('-output_fasta', required=True,
                            type=str, help='output fasta', default="none")

    parser_c5 = subparsers.add_parser("keep_best_blast_hits",
                                      help='from filterd mapping results, only keep best hit by every seqid ')
    
    parser_c5.add_argument('-input_mapping_file', required=True,
                        type=str, help='input species mapping file path (this mapping file after score filtering)', default="none")
    parser_c5.add_argument('-output_path', required=True,
                            type=str, help='output path', default="none")

    args = parser.parse_args()
    if args.subcommand == 'add_gene':
        input_gcf = args.input_gcf
        input_score = args.input_score
        temp_folder = args.temp_folder
        output_csv = args.output_csv
        time_start_s = time.time()
        add_gene(input_gcf, input_score, temp_folder, output_csv)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    elif args.subcommand == 'find_unique_count':
        sorted_score_file = args.input_sorted_score_file
        whole_umi_fastq = args.input_umi_file
        input_keys = args.input_keys
        output_unique = args.output_unique
        time_start_s = time.time()
        find_unique_count(sorted_score_file, whole_umi_fastq, input_keys, output_unique)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    elif args.subcommand == 'find_unique_count_mirna':
        sorted_score_file = args.input_sorted_score_file
        whole_umi_fastq = args.input_umi_file
        output_unique = args.output_unique
        time_start_s = time.time()
        find_unique_count_mirna(sorted_score_file, whole_umi_fastq, output_unique)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    elif args.subcommand == 'extract_sequences_from_genome':
        input_blast_result = args.input_blast_result
        input_genome = args.input_genome
        output_fasta = args.output_fasta
        time_start_s = time.time()
        extract_sequences_from_genome(input_blast_result, input_genome, output_fasta)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    elif args.subcommand == 'keep_best_blast_hits':
        input_mapping_file = args.input_mapping_file
        output_path = args.output_path
        time_start_s = time.time()
        keep_best_blast_hits(input_mapping_file, output_path)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    else:
        print("Wrong input. Check parameters")


if __name__ == "__main__":
    main()