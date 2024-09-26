import pandas as pd
import os
from os import listdir, path, makedirs
import argparse
from Bio import SeqIO

def keep_best_blast_hits(original_mapping_file, output_path):
    df_unique_mapping = pd.DataFrame()
    df_unique_mapping = find_best_blast_hits(original_mapping_file)
    df_unique_mapping.to_csv(output_path, index=False)

def find_best_blast_hits(original_mapping_file):
    df_mapping = pd.read_csv(original_mapping_file, header=None)
    df_mapping.columns = ['qseqid','sacc','sstart','send','evalue','bitscore','qcovhsp','pident']
    print(df_mapping)
    # only keep one line of highest mapping result for each qseqid
    #### Sort the DataFrame by 'qseqid' by ascending order, 'pident' first then 'qcovhsp' in descending order
    df_sorted = df_mapping.sort_values(by=['qseqid', 'pident', 'qcovhsp'], ascending=[True, False, False])
    print(f"df_sorted: {df_sorted}")
    #### Drop duplicates based on 'qseqid', keeping the first occurrence
    df_unique_mapping = df_sorted.drop_duplicates(subset='qseqid', keep='first').reset_index(drop=True)
    print(f"df_unique_mapping: {df_unique_mapping}")
    return df_unique_mapping

def filtering_file(input_mapping_file, input_fasta, output_file):
    ids = []
    with open(input_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ids.append(record.id)
    # Create a DataFrame from the list of IDs
    df_fasta = pd.DataFrame([int(x) for x in ids], columns=['qseqid'])
    print(df_fasta)

    df_unique_mapping = pd.DataFrame()
    df_unique_mapping = find_best_blast_hits(input_mapping_file)

    ### use merge inner to filter mapping file
    result = pd.merge(df_fasta, df_unique_mapping, on='qseqid', how='inner')
    print(result)

    result.to_csv(output_file, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description="get top several number of species based on ranking the mapping ratesimulate a bed file")
    parser.add_argument('-input_mapping_file', required=True,
                        type=str, help='input species mapping file path (this mapping file after score filtering)', default="none")
    parser.add_argument('-input_filtered_clean_fasta', required=True,
                        type=str, help='input clean fasta file path (this file after reads length filtering)', default="none")
    parser.add_argument('-output_file', required=True,
                        type=str, help='output resutls folder path', default="none")

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_mapping_file and args.input_filtered_clean_fasta and args.output_file:
        # time_start_s = time.time()
        filtering_file(args.input_mapping_file, args.input_filtered_clean_fasta, args.output_file)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()
