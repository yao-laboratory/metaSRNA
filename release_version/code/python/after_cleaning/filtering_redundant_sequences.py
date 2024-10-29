from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import os
import pandas as pd
import argparse
import sys
import datetime
from filter_mapping_result import find_best_blast_hits

# conda_env = os.getenv('CONDA_DEFAULT_ENV')
# if conda_env:
#     print(f"Active Conda environment: {conda_env}")
# else:
#     print("No Conda environment is active.")
def check_id_duplicate(df, column_name):
    if df[column_name].duplicated().any():
        raise ValueError(f"dataframe: '{df}' , its' column: '{column_name}' has duplicate values")

def create_umi_dataframe(umi_reads):
    data_umi = []
    for index, record in enumerate(SeqIO.parse(umi_reads, "fasta")):
        data_umi.append({'qseqid': int(record.id), 'umi': str(record.seq)})
    df_umi = pd.DataFrame(data_umi, columns=['qseqid','umi'])
    check_id_duplicate(df_umi,"qseqid")
    return df_umi

def group_sequences(df, umi_reads):
    grouped_df = df.groupby('sequence').agg(
        representative_id=('qseqid', lambda x: x.iloc[0]),
        qseqid_count=('qseqid', 'count'),
        same_seq_ids=('qseqid', lambda x: '|'.join(map(str, x))),
        **({'umi_count': ('umi', lambda x: len(set(x)))} if umi_reads not in {None, 'none', 'None'} else {})
    ).reset_index()
    return grouped_df

def filtering_fasta(seq_reads, umi_reads, csv_file, output_fasta, output_csv):
    ##first step: store fasta to a dataframe
    data_seq = []
    for index, record in enumerate(SeqIO.parse(seq_reads, "fasta")):
        data_seq.append({'qseqid': int(record.id), 'sequence': str(record.seq)})
    df_fasta = pd.DataFrame(data_seq, columns=['qseqid','sequence'])
    check_id_duplicate(df_fasta,"qseqid")
    print(df_fasta)
    name_list =[]
    ###2nd step store umi fasta to a dataframe
    if umi_reads not in {None, 'none', 'None'}:
        df_umi = create_umi_dataframe(umi_reads)
        merged_df = pd.merge(df_fasta, df_umi, on='qseqid', how='left')
        check_id_duplicate(merged_df, "qseqid")
    else:
        merged_df = df_fasta
    print("umi_reads:",umi_reads)
    print("merged_df:", merged_df)
    csv_best_hits = find_best_blast_hits(csv_file)
    print("csv_best_hits:",csv_best_hits)
    blast_result = csv_best_hits[['qseqid', 'qcovhsp', 'pident']]
    merged_df_more_information = pd.merge(merged_df, blast_result, on='qseqid', how='left')
    print("merged_df_more_information:",merged_df_more_information)
    merged_df_ordered = merged_df_more_information.sort_values(by=['qcovhsp', 'pident'], ascending=[False, False])
    merged_df_ordered.reset_index(drop=True, inplace=True)
    print("merged_df_ordered:",merged_df_ordered)
    check_id_duplicate(merged_df_ordered, "qseqid")
    grouped_df = group_sequences(merged_df_ordered, umi_reads)
    grouped_sorted_df = grouped_df.sort_values(by="representative_id")
    print("grouped_sorted_df", grouped_sorted_df)
    check_id_duplicate(grouped_df, "sequence")
    grouped_sorted_df.to_csv(output_csv, index=False, header=True)
 
    ###5th store the final fasta file
    # Create a list of SeqRecord objects
    seq_records = []
    for _, row in grouped_sorted_df.iterrows():
        seq_record = SeqRecord(Seq(row['sequence']), id=str(row['representative_id']), description="")
        seq_records.append(seq_record)
    # Write the list of SeqRecord objects to a FASTA file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")


def parse_arguments():
    parser = argparse.ArgumentParser(description=" prepare mirdeep2 reads file")
    parser.add_argument('-seq_reads', required=True,
                        type=str, help='after clean fasta file', default="none")
    parser.add_argument('-umi_reads', required=True,
                        type=str, help='umi file (whole file)', default="none")
    parser.add_argument('-csv_file', required=True,
                        type=str, help='blast_score_filter.csv in map genome folder', default="none")
    parser.add_argument('-output_fasta', required=True,
                        type=str, help='output fasta file path', default="none")
    parser.add_argument('-output_csv', required=True,
                        type=str, help='output csv file path', default="none")
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.seq_reads and args.output_fasta and args.csv_file and args.output_csv:
        # time_start_s = time.time()
        filtering_fasta(args.seq_reads, args.umi_reads, args.csv_file, args.output_fasta, args.output_csv)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()