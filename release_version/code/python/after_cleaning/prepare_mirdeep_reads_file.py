import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import os

# conda_env = os.getenv('CONDA_DEFAULT_ENV')
# if conda_env:
#     print(f"Active Conda environment: {conda_env}")
# else:
#     print("No Conda environment is active.")

def get_mirdeep_reads(input_reads, csv_file, output_file):
    print("start")
    #step1. read csv file
    df_csv = pd.read_csv(csv_file, sep=",")
    df_subset_csv = df_csv[['qseqid']]
    # print(df_subset_csv)
    print("step1 finished")
    #step2. read original reads file
    df_fasta = pd.DataFrame(columns=['qseqid','sequence'])
    for index, record in enumerate(SeqIO.parse(input_reads, "fasta")):
        df_fasta.at[index, 'qseqid'] = int(record.id)
        # df.at[index, 'id'] = record.description.rsplit("_",1)[1]
        df_fasta.at[index, 'sequence'] = str(record.seq)
    # print(df_fasta)
    df_subset_fasta = pd.merge(df_subset_csv, df_fasta, on='qseqid', how='left')
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


def main(input_reads, csv_file, output_file):
    get_mirdeep_reads(input_reads, csv_file, output_file)

if __name__ == "__main__":
    input_reads = os.getenv('PARAM1', 'default_value2')
    csv_file = os.getenv('PARAM2', 'default_value3')
    output_file = os.getenv('PARAM3', 'default_value4')
    main(input_reads, csv_file, output_file)
    print("finished producing hairpin information csv file.")