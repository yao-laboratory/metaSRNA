import re
from os import listdir, path, makedirs
import pandas as pd
import os
import numpy as np
from Bio import SeqIO

EXTRA_LEFT_GAP = 40
EXTRA_RIGHT_GAP = 40

# def find_original_seq(fasta_path, df):
#     for index, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
#         df.at[index, 'seq_id'] = record.id
#         df.at[index, 'id'] = record.description.rsplit("_",1)[1]
#         df.at[index, 'sequence'] = str(record.seq)

def extract_hairpin_loops(input_linearfold_result, input_reads, csv_file, output_file):
    df = pd.DataFrame(columns=['qseqid', 'id',
                    'start', 'end', 'dis', 'structure', 'score'])
    with open(input_linearfold_result, 'r') as file:
        file_contents = file.read()
    # print(file_contents)
    pattern = re.compile(r'>(\d+) index_(\d+)\n(.*?)(?=>|$)', re.DOTALL)

    for index, match in enumerate(pattern.finditer(file_contents)):
        # print(f'working on match: {match.group(0)}...')
        df.at[index, 'qseqid'] = int(match.group(1))
        df.at[index, 'id'] = int(match.group(2))

    # matches = pattern.findall(file_contents)
    # # print(matches)

    # #step1. store all the seqid and index in hairpin predict file to df
    # for index, match in enumerate(matches):
    #     df.at[index, 'qseqid'] = int(match[0])
    #     df.at[index, 'id'] = int(match[1])
    print("step1 finished")
    #step2. read csv file
    df_csv = pd.read_csv(csv_file, sep=",")
    #Calculate the fasta length column
    df_csv['length'] = abs(df_csv['send'] - df_csv['sstart']) + 1
    df_subset_csv = df_csv[['qseqid', 'length']]
    # print(df_subset_csv)
    df_add_len = pd.merge(df, df_subset_csv, on='qseqid', how='left')
    # print(df_add_len)
    print("step2 finished")
    #step3. read original reads file
    df_fasta = pd.DataFrame(columns=['qseqid','sequence'])
    for index, record in enumerate(SeqIO.parse(input_reads, "fasta")):
        df_fasta.at[index, 'qseqid'] = int(record.id)
        # df.at[index, 'id'] = record.description.rsplit("_",1)[1]
        df_fasta.at[index, 'sequence'] = str(record.seq)
    # print(df_fasta)
    df_add_len_and_seq = pd.merge(df_add_len, df_fasta, on='qseqid', how='left')
    # print(df_add_len_and_seq)
    print("step3 finished")
    #step4. add other data in predict haripin file which needs len information
    flag = False
    hairpin_pattern = re.compile(r'Hairpin loop.*?(\d+),\s*(\d+)\).*?:\s*([\d.]+)', re.DOTALL)
    for index, match in enumerate(pattern.finditer(file_contents)):
        # # protection code: make sure the sequence in original fasta is the sequence in the predict file
        # if index + 1 > len(df):
        #     print(f"No match found in original fasta file")
        #     break
        # # if df_add_len_and_seq.at[index, 'qseqid']  and  df_add_len_and_seq.at[index, 'id']:
        # if int(match[0]) == df_add_len_and_seq.at[index, 'qseqid'] and int(match[1]) == df_add_len_and_seq.at[index, 'id']:
        # # # protection code end
        # print(match[0])
        target_length = df_add_len_and_seq.at[index, 'length']
        # df.at[index, 'qseqid'] = int(match[0])
        # df.at[index, 'id'] = int(match[1])
        verbose_section = match[2]
        # Check if 'Hairpin loop' is in the verbose section
        if 'Hairpin loop' in verbose_section:
            # first find hairpin information (find start end value)
            hairpin_match = hairpin_pattern.findall(verbose_section)
            # print(hairpin_match)
            if hairpin_match:
                # for match in hairpin_match:
                df_add_len_and_seq.at[index, 'start'] = '|'.join(
                    str(match[0]) for match in hairpin_match)
                df_add_len_and_seq.at[index, 'end'] = '|'.join(
                    str(match[1]) for match in hairpin_match)
                df_add_len_and_seq.at[index, 'dis'] = '|'.join(str((lambda start, end: -(EXTRA_LEFT_GAP + 1 - end) if end < (EXTRA_LEFT_GAP + 1)
                                                    else start - (EXTRA_LEFT_GAP + target_length) if start > (EXTRA_LEFT_GAP + target_length)
                                                    else int(0))(int(match[0]), int(match[1]))) for match in hairpin_match)
            # else:
            #     df_add_len_and_seq.at[index, 'start'] = np.nan
            #     df_add_len_and_seq.at[index, 'end'] = np.nan
            #     df_add_len_and_seq.at[index, 'dis'] = np.nan
        # next find structure information(last line in verbose)
        lines = verbose_section.strip().split('\n')
        # print(lines)
        structure = lines[-1]
        # check if the last line in verbose is stucture or not
        if "." in structure:
            # Split the string from the right on the first occurrence of '('
            parts = structure.rsplit('(', 1)
            df_add_len_and_seq.at[index, 'structure'] = parts[0].strip()
            # print(parts[1])
            number_str = parts[1].rstrip(')').strip()
            df_add_len_and_seq.at[index, 'score'] = float(number_str)
        if flag == False:
            print("step 4 start")
        df_add_len_and_seq.iloc[[index]].to_csv(output_file, index=False, header=not flag, mode='a')
        #add header for the final file
        flag = True
        # else:
        #     print(f"No match found at index : {df_add_len_and_seq.at[index, 'id']}")
   

def main(input_linearfold_result, input_reads, csv_file, output_file):
    # Check if the output file exists, delete first
    if os.path.exists(output_file):
        os.remove(output_file)
        print(f"{output_file} has been deleted.")
    else:
        print(f"{output_file} does not exist.")
    extract_hairpin_loops(input_linearfold_result, input_reads, csv_file, output_file)

if __name__ == "__main__":
    input_linearfold_result = os.getenv('PARAM1', 'default_value1')
    input_reads = os.getenv('PARAM2', 'default_value2')
    csv_file = os.getenv('PARAM3', 'default_value3')
    output_file = os.getenv('PARAM4', 'default_value4')
    main(input_linearfold_result, input_reads, csv_file, output_file)
    print("finished producing hairpin information csv file.")