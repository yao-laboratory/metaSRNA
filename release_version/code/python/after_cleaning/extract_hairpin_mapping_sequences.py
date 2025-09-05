from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def check_id_duplicate(df, column_name):
    if df[column_name].duplicated().any():
        raise ValueError(f"dataframe: '{df}' , its' column: '{column_to_check}' has duplicate values")

def extract_seuquences(input_fasta, input_mapping, output_path):
    df_names = ["input_fasta", "input_mapping"]
    df = {name: pd.DataFrame() for name in df_names}

    fasta_data = [{'qseqid': int(record.id), 'sequence': str(record.seq)} for record in SeqIO.parse(input_fasta, "fasta")]
    df["input_fasta"] = pd.DataFrame(fasta_data, columns=['qseqid', 'sequence'])
    print("df[input_fasta]:\n",df["input_fasta"] )
    check_id_duplicate(df["input_fasta"], "qseqid")

    mirna_data = pd.read_csv(input_mapping, sep='\t', header=None, names=[
        "qseqid", "sseqid", "stitle", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bits"
    ])
    print("mirna_data['qseqid']:\n", mirna_data["qseqid"])
    df["input_mapping"] = mirna_data["qseqid"].drop_duplicates().reset_index(drop=True).to_frame()
    print("df['input_mapping']:\n", df["input_mapping"])


    merged_data = df["input_mapping"].merge(df["input_fasta"], on="qseqid", how="left")
    print(merged_data)
    final_data = merged_data["sequence"].drop_duplicates().reset_index(drop=True).to_frame()
    final_data.columns =["sequence"]
    print(final_data)
    # Write DataFrame to an Excel file
    final_data.to_csv(output_path, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract hairpin mappings sequences.")
    parser.add_argument('-input_fasta', required=True,
                        type=str, help='input fasta path', default="none")
    parser.add_argument('-input_mapping', required=True,
                        type=str, help='input mapping path', default="none")
    parser.add_argument('-output_path', required=True,
                        type=str, help='output path', default="none")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_fasta and args.input_mapping and args.output_path:
        # time_start_s = time.time()
        extract_seuquences(args.input_fasta, args.input_mapping, args.output_path)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()
    print("finished extracting hairpin mappings sequences.")

