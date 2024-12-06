from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import subprocess
import argparse
import sys
from Bio import Entrez


def union_files(folder_path, top_cnt):
    files = os.listdir(folder_path)
    # file = [file for file in files if file.endswith(
    #     '.txt')]
    # # for file in txt_files:
    file_path = os.path.join(folder_path, "blastn_refprok_evalue10.txt")
    df = []
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
    except pd.errors.ParserError:
        df = pd.read_csv(file_path, sep='\t',
                            on_bad_lines='skip', header=None)
    #1st,find top cnt higher percentage mapping staxids
    new_df = df.iloc[:, [0, -3]]  # 0 or -5 qacc=0, qseqid=-5
    new_df.columns = ['qseqid', 'staxids']
    # print(new_df)
    # Remove duplicate rows
    df_no_duplicates = new_df.drop_duplicates()
    # wanna know for each staxids, how many different qseqid number
    result = df_no_duplicates.groupby('staxids')['qseqid'].nunique()
    # print(result)
    # Sort the result from highest to lowest
    result_sorted = result.sort_values(ascending=False)
    # print(result_sorted)
    ##series to dataframe
    result_sorted = result_sorted.reset_index()
    total_df = result_sorted['staxids'].drop_duplicates().head(top_cnt)
    # print(total_df)
    ##series to dataframe again
    total_df.columns=['staxids']
    total_df = total_df.reset_index()
    #print(total_df)

    #2nd find corresponding sacc number
    compare_df = df.iloc[:, [-3, 1]]
    compare_df.columns = [ 'staxids', 'sacc']
    #print(compare_df)

    total_df_two_columns = total_df.merge(
        compare_df, on='staxids', how='inner')
    #print(total_df_two_columns)
    df_unique = total_df_two_columns.drop_duplicates(subset=['index'])
    df_unique = df_unique.drop('index', axis=1).reset_index(drop=True)
    #print(df_unique)
    return df_unique
    # end mapping file

def union_top_species_old(folder_path, top_cnt, output_folder):
    df_read = union_files(folder_path, top_cnt)
    # print(df_read)
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
    df_autoget.to_csv(path.join(output_folder, 'mapping.csv'), sep="\t", index=False)

def fetch_assembly_accession_from_dblink(sacc):
    try:
        # Fetch the sequence summary using Entrez.efetch
        handle = Entrez.efetch(db="nucleotide", id=sacc, rettype="gb", retmode="text")
        records = handle.read()  # Read the raw GenBank text
        handle.close()
        
        # Search for the Assembly accession in the DBLINK section
        for line in records.splitlines():
            if line.strip().startswith("Assembly:"):
                # Extract and return the Assembly accession
                return line.split("Assembly:")[1].strip()
        return None  # Return None if no Assembly accession is found
    except Exception as e:
        print(f"Error fetching AssemblyAccession for {sacc}: {e}")
        return None

def union_top_species(folder_path, top_cnt, output_folder):
    df_read = union_files(folder_path, top_cnt)
    # print(df_read)
    # ####using command to add corresponding gcf number
    # Get the last column of the DataFrame
    last_column = df_read[df_read.columns[1]].apply(lambda x: str(x) + ".1")
    # last_column = df_read[df_read.columns[1]].apply(lambda x: str(x))
    df_autoget = pd.DataFrame(columns=['sacc', 'gcf'])
    # Process each value in the last column
    default_gcf = ""
    for value in last_column:
        # Run the shell command to get the RefSeq for each value
        # command1 = f"module load entrez-direct/16.2"
        # command2 = f"esearch -db assembly -query {value} | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"
        # combined_command = f"{command1} && {command2}"
        # output = subprocess.run(combined_command, shell=True, capture_output=True, text=True)
        # # print(output.stdout.strip())
        # # Split the line into two values
        # values = output.stdout.strip().split(' ')
        # print(values)
        assembly_accession = fetch_assembly_accession_from_dblink(value)
        df_autoget.loc[len(df_autoget)] = {'sacc': value, 'gcf': assembly_accession if assembly_accession else default_gcf}
    print(df_autoget)
    # print(f"{value}\t{output.stdout.strip()}")
    df_autoget.to_csv(path.join(output_folder, 'mapping.csv'), sep="\t", index=False)

def parse_arguments():
    parser = argparse.ArgumentParser(description="get top several number of species based on ranking the mapping ratesimulate a bed file")
    parser.add_argument('-input_folder', required=True,
                        type=str, help='input folder name', default="none")
    parser.add_argument('-top_count', required=True,
                        type=str, help='Retain the top <n> mapping species sacc and gcf values.', default="none")
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output resutls folder path', default="none")

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_folder and args.top_count and args.output_folder:
            # time_start_s = time.time()
        union_top_species(args.input_folder, int(args.top_count), args.output_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()
