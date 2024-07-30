from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Alignment
from openpyxl.styles import Border, Side
from Bio import Entrez
import subprocess
import argparse
import sys


def union_files(folder_path):
    files = os.listdir(folder_path)
    txt_files = [file for file in files if file.endswith(
        '.txt')]
    dataframes = []
    total_df = pd.DataFrame({'staxids': []})
    compare_df_final = pd.DataFrame({'sacc': [], 'staxids': []})

    mergenew_df = pd.DataFrame({'staxids': []})
    # with pd.ExcelWriter(f) as writer:
    for file in txt_files:
        file_path = os.path.join(folder_path, file)
        #print(file_path)
        df = []
        try:
            df = pd.read_csv(file_path, sep='\t', header=None)
        except pd.errors.ParserError:
            df = pd.read_csv(file_path, sep='\t',
                             on_bad_lines='skip', header=None)

        new_df = df.iloc[:, [0, -3]]  # 0 or -5 qacc=0, qseqid=-5
        new_df.columns = ['qseqid', 'staxids']
        #print(new_df)
        # Remove duplicate rows
        df_no_duplicates = new_df.drop_duplicates()
        # wanna know for each staxids, how many different qseqid number
        result = df_no_duplicates.groupby('staxids')['qseqid'].nunique()
        # Sort the result from highest to lowest
        result_sorted = result.sort_values(ascending=False).head(20)
        #print(result_sorted)
        total_df = total_df.merge(result_sorted.reset_index()[
                                  'staxids'], on='staxids', how='outer')
        #print("final", total_df)

        compare_df = df.iloc[:, [1, -3]]
        compare_df.columns = ['sacc', 'staxids']
        compare_df_nodu = compare_df.drop_duplicates()
        compare_df_final = pd.concat(
            [compare_df_nodu, compare_df_final], axis=0)

    #print(total_df)
    # no_duplicates_total_df = total_df.drop_duplicates()
    # print(no_duplicates_total_df)
    total_df_two_columns = total_df.merge(
        compare_df_final, on='staxids', how='inner')
    #print(total_df_two_columns)
    # total_df_two_columns.drop_duplicates(subset='staxids')
    df = total_df_two_columns.drop_duplicates(subset='staxids')
    # .to_csv(path.join(folder_path, 'mapping.txt'), sep='\t', index=False, header=False)
    # df = pd.read_csv(path.join(folder_path, 'mapping.txt'),sep='\t')
    df.columns = ['staxids', 'sacc']
    #print(df)
    sorted_df = df.sort_values(by='staxids', ascending=True)
    return sorted_df
    # sorted_df.to_csv(path.join(folder_path, 'mapping.csv'),
    #                  sep="\t", index=False)
    # end mapping file


def union_top_species(folder_path):
    df_read = union_files(folder_path)
    print(df_read)
    # f_read = path.join(folder_path, 'mapping.csv')
    # f2_read = path.join(path.dirname(__file__), "../../results/blast_refprok_output/6samples_evalue10/mapping.txt")
    ###using command to add corresponding gcf number based on former mapping txt
    # df_read = pd.read_csv(f_read,sep="\t")
    # df2_read = pd.read_csv(f2_read,sep="\t")
    # # print(df_read)
    # merged = pd.merge(df_read, df2_read, on='staxids', how='outer')
    # print(merged)
    # # # Create a new column in file1 to store the corresponding values from file2
    # df['gcf_number'] = merged['gcfnumber'].fillna('null')
    # # Save the updated file1
    # df.to_csv(path.join(folder_path, 'mapping.txt'), sep="\t", index=False)


    # ####using command to add corresponding gcf number
    # Get the last column of the DataFrame
    last_column = df_read[df_read.columns[1]].apply(lambda x: str(x) + ".1")
    # last_column = df_read[df_read.columns[1]].apply(lambda x: str(x))
    df_autoget = pd.DataFrame(columns=['sacc', 'gcf'])
    # Process each value in the last column
    default_gcf = ""
    for value in last_column:
        # Run the shell command to get the RefSeq for each value
        command1 = f"module load entrez-direct/16.2"
        command2 = f"esearch -db assembly -query {value} | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"
        combined_command = f"{command1} && {command2}"
        output = subprocess.run(combined_command, shell=True, capture_output=True, text=True)
        # print(output.stdout.strip())
        # Split the line into two values
        values = output.stdout.strip().split(' ')
        # print(values)
        df_autoget.loc[len(df_autoget)] = {'sacc': value, 'gcf': values[0] if len(values) > 0 else default_gcf}
    print(df_autoget)
    # print(f"{value}\t{output.stdout.strip()}")
    df_autoget.to_csv(path.join(folder_path, 'mapping.csv'), sep="\t", index=False)

def parse_arguments():
    parser = argparse.ArgumentParser(description="get top several number of species based on ranking the mapping ratesimulate a bed file")
    parser.add_argument('-input_folder', required=True,
                        type=str, help='input folder name', default="none")

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_folder:
        # time_start_s = time.time()
        union_top_species(args.input_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()
