import pandas as pd
import os
from os import listdir, path, makedirs
import argparse
from Bio import SeqIO


def filtering_file(input_mapping_file, input_fasta, output_file):
    # folder_path = path.abspath(path.join(path.dirname(__file__), f"../../../mytuber/results/{sample_id}/blast_result/"))
    # f = path.join(folder_path, "blast_score_filter.txt'")
    # f_new = path.join(output_folder, "blast_score_final_filter.txt")
    # Open and read the FASTA file
    ids = []
    with open(input_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ids.append(record.id)
    # Create a DataFrame from the list of IDs
    df_fasta = pd.DataFrame(ids, columns=['qseqid'])
    df_fasta = df_fasta.drop_duplicates().reset_index(drop=True)
    print(df_fasta)
    df_mapping = pd.read_csv(input_mapping_file, header=None)
    df_mapping.columns = ['qseqid','sacc','sstart','send','evalue','bitscore','qcovhsp','pident']
    print(df_mapping)
    ###covert them to integer
    df_fasta['qseqid'] = pd.to_numeric(df_fasta['qseqid'], errors='coerce')  # 'coerce' will set invalid parsing as NaN
    df_mapping['qseqid'] = pd.to_numeric(df_mapping['qseqid'], errors='coerce')

    ### use merge inner to filter mapping file
    result = pd.merge(df_fasta, df_mapping, on='qseqid', how='inner')
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
