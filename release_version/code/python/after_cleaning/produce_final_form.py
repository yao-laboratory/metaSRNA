from io import StringIO
import glob
import os
from os import listdir, path, makedirs
import re
import pandas as pd
import subprocess
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def check_id_duplicate(df, column_name):
    if df[column_name].duplicated().any():
        raise ValueError(f"dataframe: '{df}' , its' column: '{column_to_check}' has duplicate values")

def produce_linearfold_final_result(df):
    # Convert 'dis' column to lists of integers
    df['dis_list'] = df['dis'].apply(lambda x: list(map(int, x.split('|'))))
    
    # keep rows any int value in 'dis_list' is between 0 and 20
    filtered_df = df[df['dis_list'].apply(lambda x: any(0 < abs(num) < 20 for num in x))]

    return filtered_df[['sequence']]

# def extract_sign(pc_str):
#     # Extract "+" or "-"
#     return pc_str.rsplit(':', 1)[-1]

# def reverse_complement(seq):
#     # Create a Seq object and return its reverse complement
#     return str(Seq(seq).reverse_complement())

# def transform_sequences(row):
#     # print("row['{flag}']\n", row[flag])
#     new_row = row.replace('U', 'T')
#     print("new_row\n", new_row)
#     # if extract_sign(row["precursor coordinate"]) == "-":
#         # Reverse complement after replacing 'U' with 'T'
#     print(reverse_complement(new_row))
#     # return new_row
    
def merge_and_match(df1, df2, new_column_name):
    merged_df = pd.merge(df1, df2, on='sequence', how='left', indicator=True)

    # 1 if ('both'), 0 otherwise
    merged_df[new_column_name] = merged_df['_merge'].apply(lambda x: 1 if x == 'both' else 0)

    # Drop the '_merge' column
    return merged_df.drop(columns=['_merge'])

# def check_substring(seq, df2_sequences):
#     return any(sub_seq in seq for sub_seq in df2_sequences)

# def merge_include(df1, df2, new_column_name):
#     # Check if sequences in df1 are in df2 and assign 1 or 0 accordingly
#     df1[new_column_name]= df1['sequence'].apply(lambda seq: int(check_substring(seq, df2['sequence'])))
#     return df1

def get_seqnames_from_mrd_file(mrd_file_path):
    identifiers = []
    pattern = re.compile(r"(seq_\d+_x\d+)")

    # Open and read the .mrd file
    with open(mrd_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            matches = pattern.findall(line)
            if matches:
                identifiers.extend(matches)

    # Display the extracted identifiers
    for seq in identifiers:
        print(seq)
    return identifiers

def get_data_from_fasta(mrd_fasta_path, seq_names_list):
    record_dict = SeqIO.to_dict(SeqIO.parse(mrd_fasta_path, "fasta"))
    # Get sequences for the given record IDs
    sequence_list = [str(record_dict[record_id].seq) for record_id in seq_names_list if record_id in record_dict]
    df = pd.DataFrame(sequence_list, columns=["sequence"])
    print("original mirdeep2 sequences:\t", df)
    df.drop_duplicates(subset=["sequence"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    print("mirdeep2 sequences after filtering redundants:\t", df)

    return df

def produce_mirdeep2_results(input_mirdeep2_folder):
    df_result = pd.DataFrame(columns=["sequence"])
    if input_mirdeep2_folder:
        # print(input_mirdeep2_folder)
        mrd_temp_folder = os.path.join(input_mirdeep2_folder,"mirdeep_runs")
        # print(mrd_temp_folder)
        mrd_file_path = glob.glob(os.path.join(mrd_temp_folder, "run_*/output.mrd"))
        print(mrd_file_path)
        if len(mrd_file_path) > 1:
            raise ValueError(f"predict_mirdeep folder has duplicate mirdeep results or no mirdeep2 result")
        elif len(mrd_file_path) == 1:
            seq_names_list = []
            seq_names_list = get_seqnames_from_mrd_file(mrd_file_path[0])
            if not seq_names_list:
                return df_result
            else:  
                mrd_fasta_path = os.path.join(input_mirdeep2_folder,"middle_results/reads_collapsed.fa")
                if os.path.exists(mrd_fasta_path):
                    print("seq_number:",seq_names_list)
                    print("seq_sequence:",get_data_from_fasta(mrd_fasta_path, seq_names_list))
                    return get_data_from_fasta(mrd_fasta_path, seq_names_list)
                else:
                    print("predict mirdeep2 tool has error, no fasta being found")
        elif len(mrd_file_path) == 0:
            print("predict mirdeep2 tool has error, so put mirdeep2 column all set to 0")
    else:
        print("do not run mirdeep2 tool, so put mirdeep2 column all set to 0")
    return df_result

def produce_form(final_fasta, input_mirdeep2_folder, linearfold_results, mirna_mapping_results, duplicates_information, output_folder):
    # Step 1: Initialize DataFrames
    df_names = ["result_df_fasta", "result_df_mirna", "result_df_linearfold", "result_df_mirdeep", "result_df_duplicates"]
    df = {name: pd.DataFrame() for name in df_names}

    # Step 2: Extract non-redundant sequences and their qseqid from final_fasta
    fasta_data = [{'qseqid': int(record.id), 'sequence': str(record.seq)} for record in SeqIO.parse(final_fasta, "fasta")]
    df_fasta= pd.DataFrame(fasta_data)
    check_id_duplicate(df_fasta, "qseqid")
    # check_id_duplicate(df_fasta, "sequence")
    df["result_df_fasta"] = df_fasta['sequence'].to_frame()
    print("df['result_df_fasta']:\n", df["result_df_fasta"])

    # Step 3: Extract unique qseqid from miRNA mapping file
    df["result_df_mirna"] = pd.read_csv(mirna_mapping_results, sep='\t')
    print("df['result_df_mirna']:\n", df["result_df_mirna"])
    
    # Step 4: Extract all information from duplicates information file
    df["result_df_duplicates"] = pd.read_csv(duplicates_information, sep=',')
    print("df['result_df_duplicates']:\n", df["result_df_duplicates"])

    # Step 5: Extract qseqid, sequence, distance from LinearFold prediction result
    if linearfold_results:
        linearfold_df = pd.read_csv(linearfold_results, sep=',')
        df["result_df_linearfold"] = produce_linearfold_final_result(linearfold_df)
        # check_id_duplicate(df["result_df_linearfold"], "sequence")
        print("df['result_df_linearfold']:\n", df["result_df_linearfold"])
    else:
        print("do not have linearfold results, so put linearfold column all set to 0")
        df["result_df_linearfold"] = pd.DataFrame(columns=["sequence"])

    # Step 6: Extract sequence from miRDeep2 prediction result
    df["result_df_mirdeep"] = produce_mirdeep2_results(input_mirdeep2_folder)
                
    #merge them together
    for df_name in df_names:
        check_id_duplicate(df[df_name], "sequence")
    merged_mirna = merge_and_match(df["result_df_fasta"], df["result_df_mirna"],"mirBase")
    merged_linearfold = merge_and_match(merged_mirna, df["result_df_linearfold"], "linearfold")
    merged_mirdeep = merge_and_match(merged_linearfold, df["result_df_mirdeep"],"mirdeep2")
    merged_final = merged_mirdeep.merge(df["result_df_duplicates"], on="sequence", how="left")
    print("merged_final:\n", merged_final)
    # Write DataFrame to an Excel file
    print(output_folder)
    output_path = os.path.join(output_folder, "final_form.csv")
    print(output_path)
    merged_final.to_csv(output_path, index=False)
        

def parse_arguments():
    parser = argparse.ArgumentParser(description="get top several number of species based on ranking the mapping ratesimulate a bed file")
    parser.add_argument('-input_final_fasta', required=False,
                        type=str, help='Path to the clean FASTA file (after filtering length and removing duplicates) ')
    parser.add_argument('-input_mirna_mapping_results', required=False,
                        type=str, help='Path to the miRNA mapping score file.')
    parser.add_argument('-duplicates_information', required=False,
                        type=str, help='Path to file about duplicate sequences information.')
    parser.add_argument('-input_mirdeep2_folder', required=False,
                        type=str, help='Folder contains mirdeep2 prediction result.')
    parser.add_argument('-input_linearfold_results', required=False,
                        type=str, help='Path to file about linearfold prediction results.')
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output resutls folder path')

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_final_fasta and args.input_mirna_mapping_results and args.duplicates_information and args.output_folder:
            # time_start_s = time.time()
        produce_form(args.input_final_fasta, args.input_mirdeep2_folder, args.input_linearfold_results, args.input_mirna_mapping_results, args.duplicates_information, args.output_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()