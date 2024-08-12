from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Alignment
from openpyxl.styles import Border, Side
import argparse
import sys


def find_species_classfication(folder_path, total_sequence_number, output_folder):
    file = path.abspath(path.join(folder_path, "blastn_refprok_evalue10.txt"))
    file_path = os.path.join(folder_path, file)
    # print(file_path)
    f_output = path.join(output_folder, 'refprok_species_classification_analysis.csv')
    df = []
    try:
        df = pd.read_csv(file_path, sep='\t', header= None)
    except pd.errors.ParserError:
        df = pd.read_csv(file_path, sep='\t', on_bad_lines='skip', header= None)

    new_df = df.iloc[:, [0, -3]]
    new_df.columns = ['qacc', 'staxids']
    # print(new_df)
    # Remove duplicate rows
    df_no_duplicates = new_df.drop_duplicates() 
    # print(df_no_duplicates)
    result = new_df.groupby('staxids')['qacc'].nunique()
    # result.columns = ['staxids', 'unique_qacc_count']
    # print(result)

    names_df = df.iloc[:, [-3,-2,-1]]
    names_df.columns = ['staxids', 'sscinames','sblastnames']
    no_duplicates_name_df = names_df.drop_duplicates() 
    #print(no_duplicates_name_df)
    # Sort the result from highest to lowest
    result_sorted = result.sort_values(ascending=False)
    # Merge the two DataFrames on the 'key' column
    merge_df = pd.merge(result_sorted, no_duplicates_name_df, on='staxids')
    merge_df.columns = ['staxids', 'unique_qacc_count', 'sscinames', 'sblastnames']
    merge_df['percentage'] = merge_df['unique_qacc_count'] / total_sequence_number * 100
    merge_df['percentage'] = merge_df['percentage'].apply(lambda x: '{:.2f}%'.format(x))
    # print(merge_df)
    # print(merge_df)
    # Count number of unique numbers in column 'A'
    # Append a new row

    unique_count_A = new_df['qacc'].nunique()/total_sequence_number * 100
    # print(unique_count_A)
    # Create a new DataFrame with one column
    new_data = {'blast_sequences_num':[new_df['qacc'].nunique()], 'total_afterclean_biggerthan_12seq_number':[total_sequence_number], 'blast_sequences_percentage': [unique_count_A ]}
    oneline_df = pd.DataFrame(new_data)
    oneline_df['blast_sequences_percentage'] = oneline_df['blast_sequences_percentage'].apply(lambda x: '{:.2f}%'.format(x))
    # Concatenate the original DataFrame with the new DataFrame
    final_df = pd.concat([merge_df, oneline_df], ignore_index=True)
    # print(final_df)
    final_df.to_csv(f_output, index=False)
    # with pd.ExcelWriter(f_output) as writer:
    #     # # Write DataFrame to an Excel file
    #     final_df.to_excel(writer, sheet_name=file, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description="get all mapping species and their information from refprok database")
    parser.add_argument('-input_folder', required=True,
                        type=str, help='input mapping resutls folder path', default="none")
    parser.add_argument('-total_sequence_number', required=True,
                        type=str, help='fasta file after clean and its'' total sequence number ', default="none")
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output resutls folder path', default="none")

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_folder and args.total_sequence_number and args.output_folder:
        # time_start_s = time.time()
        find_species_classfication(args.input_folder, int(args.total_sequence_number), args.output_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()