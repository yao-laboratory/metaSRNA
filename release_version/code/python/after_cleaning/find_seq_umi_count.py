import re
from io import StringIO
import os
import pandas as pd
import numpy as np
import time
import csv
import math


# whole_umi_fastq = "final_whole_umi.fastq"
# sorted_score_file = "blast_score_filter_add_gene.csv"
# output_unique_csv = "output_unique_seq_umi.csv"

def calculate_uniseq(input_string, sorted_score_file):
    tags = input_string.split()
    #print(tags)
    
    #start copy 7 columns to a new dataframe
    sorted_score_df = pd.read_csv(sorted_score_file, sep=",",header = None, dtype={8: str})
    #inital the sorted_score_new which has ID, KEY1,KEY2...
    sorted_score_new = pd.DataFrame()
    sorted_score_new['id'] = sorted_score_df.iloc[:, 0].copy()
    
    
    #set keys become new columns, initiate the value to nans
    for i in range(len(tags)):
        sorted_score_new[tags[i]] = np.nan
        sorted_score_new[tags[i]] = sorted_score_new[tags[i]].astype(str)
    #print(sorted_score_new)
    
    
    
    # seperate the key and value from the original dataframe, then put key-values to the different columns to the new dataframe
    #gene_id "ECDH1ME8569_RS00465"; transcript_id ""; gbkey "Gene"; gene "murC"; gene_biotype "protein_coding"; locus_tag "ECDH1ME8569_RS00465"; old_locus_tag "ECDH1ME8569_0088"
    pattern = r'(\w+)\s"([^"]*)";' 
    
    for index, row in sorted_score_df.iterrows():
        #find each row's last string with bio information
        if isinstance(row.iloc[-1], str):
            #print(sorted_score.iloc[-1])
            #match = re.search(pattern, sorted_score[-1])
            value = re.findall(pattern, row.iloc[-1])
            # print(value)
            result = dict(value)
            # print(result)
            # Check if a match is found
            # only calculate the first row of same sequence number 
            if result:
                for i in range(len(tags)):
                    sorted_score_new.at[index, tags[i]] = result.get(tags[i])
#             else:
#                 print(index)
#                 print('No match found')
    # sorted_score_new = sorted_score_new.fillna(0)
    print(sorted_score_new)
    
    #start create dataframe seqnumlist for id, seq_count and seq_list 
    group_list = []
    for i in range(len(tags)):
        group_list.append(tags[i])
        
    seq_count_column = sorted_score_new.drop_duplicates().groupby(group_list).count()
    seqnumlist = pd.DataFrame()
    seqnumlist= seq_count_column
    seqnumlist.columns = ['unique_seq_count']
    seqnumlist['seq_num_list'] = sorted_score_new.drop_duplicates().groupby(group_list)['id'].apply(list)
    print(seqnumlist)
    return seqnumlist

def calculate_uniseq_mirna(sorted_score_file):
    df = pd.read_csv(sorted_score_file, sep="\t",header = None, index_col=False)
    mirna_df = df.iloc[:, [0,2]]
    mirna_df.columns=['id','gene_biotype']
    seq_count_column = mirna_df.drop_duplicates().groupby('gene_biotype').count()
    print(seq_count_column)
    seqnumlist = pd.DataFrame()
    seqnumlist= seq_count_column
    seqnumlist.columns = ['unique_seq_count']
    seqnumlist['seq_num_list'] = mirna_df.drop_duplicates().groupby(mirna_df['gene_biotype'])['id'].apply(list)
    print(seqnumlist)
    return seqnumlist

#make dictionary for seq_id as key, umi sequence as dictionary
def calculate_seq_to_uniumi(whole_umi_fastq):
    #start calculate unique umi count:
    # first step: create umi dictionary
    whole_umi_dict = {}
    if whole_umi_fastq not in {None, 'none', 'None'}:
        with open(whole_umi_fastq, 'r') as input_umi_file:
            while True:
                lines = []
                for i in range(2):
                    line = input_umi_file.readline()
                    if (not line):
                        break
                    lines.append(line.strip())
                if (not lines or len(lines) <= 1):
                    input_umi_file.close()
                    break
                if len(lines) == 2 and lines[0].find('@') != -1:
                    filtered_string = lines[0].replace('@', '')
                    seq_num = int(filtered_string)
                    #print(seq_num)
                    if seq_num not in whole_umi_dict:
                        whole_umi_dict[seq_num] = ""
                    whole_umi_dict[seq_num] = lines[1]            
    return whole_umi_dict

#add unique_umi_count to dataframe seqnumlist then save to the csv file
def save_seq_umi(whole_umi_dict, output_unique_csv, seqnumlist):
    # Write the modified data back to the CSV file with the new column
    # print(whole_umi_dict)
    # print(genetype_seqnumlist_dict)
    if whole_umi_dict:
        for index, row in seqnumlist.iterrows():
            # print(row['seq_num_list'])
            unique_seqnum_set = sorted(row['seq_num_list'])
            # print("bbbb",unique_seqnum_set)
            umi_list = []
            for unique_seq_num in unique_seqnum_set:
                # print("ccc",unique_seq_num)
                umi_list.append(whole_umi_dict[int(unique_seq_num)])
            # print(umi_list)
            unique_umi_list = set(umi_list)
            seqnumlist.at[index, 'unique_umi_count'] = len(unique_umi_list)
            
    print(seqnumlist)
    df_filtered = pd.DataFrame()
    df_filtered = seqnumlist.drop('seq_num_list', axis=1)
    df_filtered.to_csv(output_unique_csv, sep=',', index = True)


def find_unique_count(sorted_score_file, whole_umi_fastq, input_keys, output_unique):
    seqnumlist = calculate_uniseq(input_keys, sorted_score_file)
    whole_umi_dict = calculate_seq_to_uniumi(whole_umi_fastq)
    save_seq_umi(whole_umi_dict, output_unique, seqnumlist)

def find_unique_count_mirna(sorted_score_file, whole_umi_fastq, output_unique):
    seqnumlist = calculate_uniseq_mirna(sorted_score_file)
    whole_umi_dict = calculate_seq_to_uniumi(whole_umi_fastq)
    save_seq_umi(whole_umi_dict, output_unique, seqnumlist)


# def main():
#     find_unique_count("blast_score_filter_add_gene.csv", "final_whole_umi.fastq", "output_unique_seq_umi.csv")
#     python3 $code_path/blast_process_main.py find_unique_count -input_sorted_score_file $blast_score/blast_score_filter_add_gene_1db_${bacteria}.csv -input_umi_file $results/run_result/final_whole_umi.fastq -input_keys "gene_biotype" -output_unique $blast_score/output_unique_biotype_1db_${bacteria}.csv
# if __name__ == "__main__":
#     main()
#     print("finished adding.")