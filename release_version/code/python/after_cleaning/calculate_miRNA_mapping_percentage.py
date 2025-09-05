from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import argparse
import sys


def calculate_percentage(total_sequence_number, input_path, output_folder):
    f = path.join(output_folder, 'hairpinrna_analysis.csv')
    total_lines = []
    df = []
    # try:
    df = pd.read_csv(input_path, sep='\t', header= None)
    # except pd.errors.ParserError:
    #     df = pd.read_csv(file_path, sep=',', on_bad_lines='skip', header= None)
    
    new_df = df.iloc[:, [0]]
    new_df.columns = ['sacc']
    # print(new_df)
    result = new_df['sacc'].nunique()
    print(new_df['sacc'].drop_duplicates())
    oneline = {'unique_sacc_number':result, 'total_sequences_number(after_clean)': total_sequence_number, 'percentage': result / total_sequence_number * 100, 'file_name':path.abspath(input_path)}
    total_lines.append(oneline)

    total_lines = pd.DataFrame(total_lines)
    total_lines['percentage'] = total_lines['percentage'].apply(lambda x: '{:.2f}%'.format(x))
    # print(total_lines.shape[1])
    # print(total_lines)
    # Write DataFrame to an Excel file
    total_lines.to_csv(f, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate the percentage of sequences mapped to the hairpin database")
    parser.add_argument('-total_seq_number', required=True,
                        type=int, help='input path', default="none")
    parser.add_argument('-input_path', required=True,
                        type=str, help='input path', default="none")
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output folder', default="none")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_path and args.output_folder:
        # time_start_s = time.time()
        calculate_percentage(args.total_seq_number, args.input_path, args.output_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()
    print("finished calculating mirna mapping percentage.")

