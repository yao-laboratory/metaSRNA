from io import StringIO
import os
from os import listdir, path, makedirs
import re
import pandas as pd
import argparse
# from Bio.Seq import Seq
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
def check_id_duplicate(df, column_name):
    if df[column_name].duplicated().any():
        raise ValueError(f"dataframe: '{df}' , its' column: '{column_name}' has duplicate values")

def build_all_sequences_results(df_result):
    expanded_rows = []
    for _, row in df_result.iterrows():
        numbers = row['same_seq_ids'].split('|')
        for number in numbers:
            new_row = row.copy()
            new_row['same_seq_ids'] = number
            expanded_rows.append(new_row)

    new_df = pd.DataFrame(expanded_rows)
    return new_df

def save_all_data_to_bedfile(all_data_form, bed_file_path):
    # store data first
    data = []
    # ("chr1", 2000, 6000, "Feature2", 0, "-"),
    for index, row in all_data_form.iterrows():
        chrome = row["qseqid"]
        start_axis = row["sstart"]
        end_axis = row["send"]
        sequence = row["sequence"]
        data.append((chrome, start_axis,
                    end_axis, sequence, "", "", "", ""))

    # write data to the bed file
    with open(bed_file_path, "w") as file:
        for line in data:
            file.write("\t".join(map(str, line)) + "\n")

def produce_form(final_form, mapping_with_genes_file, output_folder):
    # Step 1: Initialize DataFrames
    df_names = ["result_form", "result_mapping", "result_merge", "result_all_sequences","result_merge_all", "result_merge_all_save"]
    df = {name: pd.DataFrame() for name in df_names}
    #step1
    df["result_form"] = pd.read_csv(final_form, sep=',')
    df["result_form"] = df["result_form"].rename(columns={'representative_id': 'qseqid'})
    print("df['result_form']:\n", df["result_form"])
    #step2
    df_mapping = pd.read_csv(mapping_with_genes_file, sep=',')
    df_mapping.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "coverage", "identity", "overlap_gene"]
    df['result_mapping'] = df_mapping[["qseqid", "coverage", "identity", "overlap_gene", "sstart", "send"]]
    print("df['result_mapping']:\n", df["result_mapping"])
    #step3
    df["result_merge"] = pd.merge(df["result_form"], df["result_mapping"], on='qseqid', how='left')
    move_column = df["result_merge"].pop("same_seq_ids")
    df["result_merge"]["same_seq_ids"] = move_column
    new_column_order = ['qseqid'] + [col for col in df["result_merge"].columns if col != 'qseqid']
    df["result_merge"]= df["result_merge"][new_column_order]
    check_id_duplicate(df['result_merge'], "qseqid")
    print("df['result_merge']:\n", df["result_merge"])
    
    output_path_1 = os.path.join(output_folder, "representative_sequence_results.csv")
    print(output_path_1)
    df["result_merge"].to_csv(output_path_1, index=False)

    #step4
    df["result_all_sequences"] = build_all_sequences_results(df["result_merge"])
    print("df['result_all_sequences'].columns:\n", df["result_all_sequences"].columns)
    df["result_all_sequences"] = df["result_all_sequences"].drop(['qseqid', 'coverage', 'identity', 'overlap_gene', 'sstart', 'send'], axis=1)
    df["result_all_sequences"] = df["result_all_sequences"].rename(columns={'same_seq_ids': 'qseqid'})
    df["result_all_sequences"]['qseqid'] = df["result_all_sequences"]['qseqid'].astype(int)
    print("df['result_all_sequences']:\n", df["result_all_sequences"])
    df["result_merge_all"] = pd.merge(df["result_all_sequences"], df["result_mapping"], on='qseqid', how='left')
    ##make qseqid to the first column
    new_column_order = ['qseqid'] + [col for col in df["result_merge_all"].columns if col != 'qseqid']
    df["result_merge_all"] = df["result_merge_all"][new_column_order]
    print("df['result_merge_all'].columns:\n", df['result_merge_all'].columns)
    print("df['result_merge_all']:\n", df["result_merge_all"])
    check_id_duplicate(df["result_merge_all"], "qseqid")

    df["result_merge_all_save"] = df["result_merge_all"].copy()
    df["result_merge_all_save"].drop(['sstart', 'send'], axis=1, inplace=True)
    print("df['result_merge_all_save']:\n", df["result_merge_all_save"])
    output_path_2 = os.path.join(output_folder, "all_sequence_results.csv")
    print(output_path_2)
    df["result_merge_all_save"].to_csv(output_path_2, index=False)

    #step5
    output_path_3 = os.path.join(output_folder, "representative_sequence_results.bed")
    print(output_path_3)
    save_all_data_to_bedfile(df["result_merge"],output_path_3)
    output_path_4 = os.path.join(output_folder, "all_sequence_results.bed")
    print(output_path_4)
    save_all_data_to_bedfile(df["result_merge_all"],output_path_4)

def parse_arguments():
    parser = argparse.ArgumentParser(description="build additional step")
    parser.add_argument('-input_final_form', required=False,
                        type=str, help='final form')
    parser.add_argument('-input_mapping_with_genes_file', required=False,
                        type=str, help='Path to the miRNA mapping score file with gene information.')
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output resutls folder path')

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_final_form and args.input_mapping_with_genes_file and args.output_folder:
        # time_start_s = time.time()
        produce_form(args.input_final_form, args.input_mapping_with_genes_file, args.output_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()