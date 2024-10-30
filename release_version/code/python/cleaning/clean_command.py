import argparse
import time
from clean_main_step import clean_main_step
from clean_main_step import after_clean_step
from clean_main_step import clean_N_seqs
    
def main():
    parser = argparse.ArgumentParser(prog='Cleanfastaq')

    subparsers = parser.add_subparsers(
        dest='subcommand', help='Sub Commands Help')
    # add sub command
    parser_cm = subparsers.add_parser("clean_fastq",
                                      help='')
    parser_cm.add_argument('-input', required=True,
                           type=str, help='input file: fastq format', default="none")
    parser_cm.add_argument('-fa_format', required=True,
                       type=str, help='format', default="none")
    parser_cm.add_argument('-output_filename', required=True,
                           type=str, help='output file main name', default="none")
    parser_cm.add_argument('-fault_tolerance', required=True,
                           type=str, help='', default="none") 
    parser_cm.add_argument('-tail_incomplete_tolerance', required=True,
                       type=str, help='', default="none") 
    
    parser_ac = subparsers.add_parser("after_clean",
                                      help='')
    parser_ac.add_argument('-input', required=True,
                           type=str, help='input file: fastq format', default="none")
    parser_ac.add_argument('-output_path', required=True,
                           type=str, help='output file path', default="none")

    parser_cns = subparsers.add_parser("clean_N_seqs",
                                      help='')
    parser_cns.add_argument('-input', required=True,
                           type=str, help='input file: fastq format', default="none")
    parser_cns.add_argument('-output_path', required=True,
                           type=str, help='output file path', default="none")

    args = parser.parse_args()
    if args.subcommand == 'clean_fastq':
        input_file = args.input
        output_filename = args.output_filename
        fa_format = args. fa_format
        tolerance = args.fault_tolerance
        tail_tolerance = args.tail_incomplete_tolerance
        time_start_s = time.time()
        clean_main_step(input_file,output_filename,fa_format, int(tolerance), int(tail_tolerance))
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    elif args.subcommand == 'after_clean':
        input_file = args.input
        output_file = args.output_path
        time_start_s = time.time()
        after_clean_step(input_file, output_file)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    elif args.subcommand == 'clean_N_seqs':
        input_file = args.input
        output_file = args.output_path
        time_start_s = time.time()
        clean_N_seqs(input_file, output_file)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    else:
        print("Wrong input. Check parameters")


if __name__ == "__main__":
    main()