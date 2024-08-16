import re
from os import listdir, path, makedirs
import pandas as pd
import os
import numpy as np
from Bio import SeqIO
import argparse
import sys

EXTRA_LEFT_GAP = 40
EXTRA_RIGHT_GAP = 40

def calculate_dis(row):
    start_values = list(map(int, row['start'].split('|')))
    end_values = list(map(int, row['end'].split('|')))
    target_length = int(row['length'])
    dis = '|'.join(str((lambda start, end: -(EXTRA_LEFT_GAP + 1 - end) if end < (EXTRA_LEFT_GAP + 1) \
                else start - (EXTRA_LEFT_GAP + target_length) if start > (EXTRA_LEFT_GAP + target_length) \
                else int(0))(start, end)) for start, end in zip(start_values, end_values))
    return dis

def extract_hairpin_loops(input_linearfold_result, input_reads, csv_file, output_file):
    #step1. read linearfold results
    with open(input_linearfold_result, 'r') as file:
        file_contents = file.read()
    # print(file_contents)
    seq_pattern = re.compile(r'>(\d+) index_(\d+)\n(.*?)(?=>|$)', re.DOTALL)
    program_data = []
    for index, match in enumerate(seq_pattern.finditer(file_contents)):
        #code protection
        try:
            if match.group(1) is None or match.group(2) is None or match.group(3) is None:
                continue
        except IndexError:
            continue
        start = end = dis = structure = score = np.nan
        verbose_section = match.group(3)
        hairpin_pattern = re.compile(r'Hairpin loop.*?(\d+),\s*(\d+)\).*?:\s*([\d.]+)', re.DOTALL)
        if 'Hairpin loop' in verbose_section:
            # first find hairpin information (find start end value)
            hairpin_match = hairpin_pattern.findall(verbose_section)
            if hairpin_match:
                starts, ends = zip(*[(str(match[0]), str(match[1])) for match in hairpin_match])
                start = '|'.join(starts)
                end = '|'.join(ends)
        lines = verbose_section.strip().split('\n')
        # print(lines)
        if lines:
            structure = lines[-1]
            # check if the last line in verbose is stucture or not
            if "." in structure:
                # Split the string from the right on the first occurrence of '('
                parts = structure.rsplit('(', 1)
                try:
                    if len(parts) != 2:
                        continue
                except IndexError:
                    continue
                structure = parts[0].strip()
                number_str = parts[1].rstrip(')').strip()
                score = float(number_str)
        program_data.append({'qseqid': int(match.group(1)),'id': int(match.group(2)), 'start': start, 'end': end, 'structure': structure, 'score': score})
    df = pd.DataFrame(program_data, columns=['qseqid', 'id',
                    'start', 'end', 'dis', 'structure', 'score'])
    # print(df)
    print("step1 finished")

    #step2. read csv file
    df_csv = pd.read_csv(csv_file, sep=",")
    df_csv.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "qcovhsp", "pident"]
    #Calculate the fasta length column
    df_csv['length'] = abs(df_csv['send'] - df_csv['sstart']) + 1
    df_subset_csv = df_csv[['qseqid', 'length']].drop_duplicates()
    # print(df_subset_csv)
    df_add_len = pd.merge(df, df_subset_csv, on='qseqid', how='left')
    # print(df_add_len)
    print("step2 finished")

    #step3. read original reads file
    data = []
    for record in SeqIO.parse(input_reads, "fasta"):
        data.append({'qseqid': int(record.id),'sequence': str(record.seq)})
    df_fasta = pd.DataFrame(data, columns=['qseqid','sequence'])
    # print(df_fasta)
    df_add_len_and_seq = pd.merge(df_add_len, df_fasta, on='qseqid', how='left')
    # print(df_add_len_and_seq)
    print("step3 finished")

    #step4. calculate distances
    #filter nan first
    final_df = df_add_len_and_seq.dropna(subset=['start', 'end', 'length'])
    # Apply the function to each row
    final_df = final_df.copy()
    final_df['dis'] = final_df.apply(calculate_dis, axis=1)
    final_df.to_csv(output_file, index=False)
    print("step 4 finished")
    # print(final_df)
    # Check for any NaN in the entire DataFrame
    has_nan = final_df.isna().any().any()
    print(f"Does the DataFrame have any NaN? {has_nan}")

def parse_arguments():
    parser = argparse.ArgumentParser(description=" produce hairpin information csv file")
    parser.add_argument('-input_linearfold_result', required=True,
                        type=str, help='linearfold output path', default="none")
    parser.add_argument('-input_reads', required=True,
                        type=str, help='after clean fasta file', default="none")
    parser.add_argument('-csv_file', required=True,
                        type=str, help='mapping species database csv path', default="none")
    parser.add_argument('-output_file', required=True,
                        type=str, help='output file path', default="none")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if os.path.exists(args.output_file):
        os.remove(args.output_file)

    if args.input_linearfold_result and args.input_reads and args.csv_file and args.output_file:
        # time_start_s = time.time()
        extract_hairpin_loops(args.input_linearfold_result, args.input_reads, args.csv_file, args.output_file)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()