import re
from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import numpy as np
import time
import csv

# import scipy
# from skbio import TreeNode, read

# gtf_file = "gene.gtf"
# score_file = "score_filter.txt"
# file_path = "blast_score_filter_add_gene.csv"

def print_using_time(index_s, start_time, contig):
    if index_s % 100000 == 0:
        end = time.time() 
        use = end - start_time
        total_line = (index_s + 1) if index_s == 0 else index_s
        print(f"contig name is {contig}, and already processed this contig's {total_line} lines in score file, total used time is {use} seconds")

def add_gene_information_by_every_contig(mapping_file_part, gtf_df_part, output_file, start_time, contig):
    gtf_index = 0
    gtf_index_previous = 0
    for index_s, row_s in mapping_file_part.iterrows():
        #print time for each 10000 rows processed
        print_using_time(index_s, start_time, contig)
        #finish one score line's gene finding
        gtf_index = finding_one_score_line_gene(index_s, row_s, mapping_file_part, gtf_df_part, gtf_index_previous)
        gtf_index_previous = gtf_index
    
    mapping_file_part.to_csv(output_file, sep=',', mode='a', header=None, index=False, quoting=csv.QUOTE_NONE, escapechar='\\')
        
def finding_one_score_line_gene(index_s, row_s, sorted_score_df, gtf_df, gtf_index_previous):
     for index_g, row_g in gtf_df[gtf_index_previous:].iterrows():
        # print(index_g,gtf_index_previous)
        # if row_g["source"] == "RefSeq" and row_g["feature"] == "gene":
        start_pos = int(row_s["sstart"]) if int(row_s["sstart"]) <= int(row_s["send"]) else int(row_s["send"])
        end_pos =  int(row_s["send"]) if int(row_s["sstart"]) <= int(row_s["send"]) else int(row_s["sstart"])

        #if score start&end axies in the gene area,then take the gene information to this score line 
        if int(row_g["start"]) <= start_pos and end_pos <= int(row_g["end"]):
            #sorted_score_df.at[index_s,8] = str(row_g[8])
            sorted_score_df.at[index_s, "gene_information"]= row_g["attribute"]
            return index_g
        # if score area on the left, gene inoformation area on the right, 
            # means before not in the area,cannot find the gene name, 
            # should return the none found result and also return the inital place for this searching 
        if  start_pos < int(row_g["start"]):
            return gtf_index_previous
    # return gtf_index_previous


def add_gene(input_gtf, input_score, temp_folder, output_csv):
    gtf_df = pd.read_csv(input_gtf, sep="\t", header=None)
    gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf_df = gtf_df[(gtf_df["source"] == "RefSeq") & (gtf_df["feature"] == "gene")]
    #gtf needs reorder
    gtf_df.sort_values(by=['seqname', 'start', 'end'], ascending=[True, True, True], inplace=True)
    gtf_df.reset_index(drop=True, inplace=True)
    gtf_df["attribute"] = gtf_df["attribute"].str.replace(',', ';')
    print(gtf_df.head())

    mapping_score_df = pd.read_csv(input_score, sep=",")
    mapping_score_df.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "qcovhsp", "pident"]
    mapping_score_df["gene_information"] = np.nan
    #mapping score file needs reorder
    mapping_score_df.sort_values(by=['sacc', 'sstart', 'send'], ascending=[True, True, True], inplace=True)
    mapping_score_df.reset_index(drop=True, inplace=True)
    print("mapping_score_df:", mapping_score_df.head())
    # temp_csv = os.path.join(temp_folder, "blast_score_filter_ordered_temp.csv")
    # temp.to_csv(temp_csv, sep=',', header=True, index=False)
    start_time = time.time() 

    temp_final_csv = os.path.join(temp_folder, "blast_score_filter_add_gene_temp.csv")
    with open(temp_final_csv, 'w') as out_file:
        gtf_groups = gtf_df.groupby('seqname')
        mapping_score_groups = mapping_score_df.groupby('sacc')
        for contig, mapping_score_group in mapping_score_groups:
            if contig in gtf_groups.groups:
                ##for each contig pairs, need to reset index
                mapping_score_group = mapping_score_group.reset_index(drop=True)
                gtf_group = gtf_groups.get_group(contig).reset_index(drop=True)
                add_gene_information_by_every_contig(mapping_score_group, gtf_group, out_file, start_time, contig)

            else:
                print(f"Contig {contig} in score file not found in gtf")

        # add_gene_information_by_every_contig(temp, gtf_df, output_csv)
    temp_final_df = pd.read_csv(temp_final_csv, sep=",", header=None)
    temp_final_df.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "qcovhsp", "pident", "gene_information"]
    temp_final_df.sort_values(by=['qseqid', 'sstart', 'send'], ascending=[True, True, True], inplace=True)
    temp_final_df.reset_index(drop=True, inplace=True)
    temp_final_df["gene_information"] = temp_final_df["gene_information"].fillna('nan')
    temp_final_df.to_csv(output_csv, sep=",", header=None, index=False, quoting=csv.QUOTE_NONE, escapechar='\\')

    # if os.path.exists(temp_csv):
        #     os.remove(temp_csv)
    if os.path.exists(temp_final_csv):
        os.remove(temp_final_csv)
        

