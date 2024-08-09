from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import os
import pandas as pd
import argparse
import sys
print(sys.executable)

# conda_env = os.getenv('CONDA_DEFAULT_ENV')
# if conda_env:
#     print(f"Active Conda environment: {conda_env}")
# else:
#     print("No Conda environment is active.")

def get_mirdeep_reads(input_reads, csv_file, output_file):
    print("start")
    #step1. read csv file
    df_csv = pd.read_csv(csv_file, sep=",")
    df_csv.columns=['qseqid', 'sacc', 'sstart', 'send', 'evalue', 'bitscore', 'qcovhsp', 'pident']
    # print(df_csv)
    df_subset_csv = df_csv['qseqid'].drop_duplicates()
    # print(df_subset_csv)
    print("step1 finished")

    #step2. read original reads file
    df_fasta = pd.DataFrame(columns=['qseqid','sequence'])
    for index, record in enumerate(SeqIO.parse(input_reads, "fasta")):
        df_fasta.at[index, 'qseqid'] = int(record.id)
        # df.at[index, 'id'] = record.description.rsplit("_",1)[1]
        df_fasta.at[index, 'sequence'] = str(record.seq)
    # print(df_fasta)
    ### use inner merge cause after cleaning fasta filtering 17 and mapping with after cleaning fasta filtering 12 have different data, needs finding overlap
    df_subset_fasta = pd.merge(df_subset_csv, df_fasta, on='qseqid', how='inner')
    # print(df_subset_fasta)
    print("step2 finished")

    #step3. store dataframe as a fasta
    # Convert DataFrame rows to SeqRecord objects, ensure sequences and id are strings
    records = [
        SeqRecord(Seq(str(row['sequence'])), id=str(row['qseqid']), description="")
        for index, row in df_subset_fasta.iterrows()
    ]

    # Write the records to the FASTA file
    SeqIO.write(records, output_file, "fasta")
    print("step3 finished")


def parse_arguments():
    parser = argparse.ArgumentParser(description=" prepare mirdeep2 reads file")
    parser.add_argument('-input_reads', required=True,
                        type=str, help='after clean fasta file', default="none")
    parser.add_argument('-csv_file', required=True,
                        type=str, help='mapping species database csv path', default="none")
    parser.add_argument('-output_file', required=True,
                        type=str, help='output file path', default="none")
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.input_reads and args.csv_file and args.output_file:
        # time_start_s = time.time()
        get_mirdeep_reads(args.input_reads, args.csv_file, args.output_file)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()