import argparse
import time
from clean_main_step import clean_main_step
    
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

    else:
        print("Wrong input. Check parameters")


if __name__ == "__main__":
    main()