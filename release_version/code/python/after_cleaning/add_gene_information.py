import re
from io import StringIO
import os
import pandas as pd
import numpy as np
import time
import csv

# import scipy
# from skbio import TreeNode, read

# gtf_file = "gene.gtf"
# score_file = "score_filter.txt"
# file_path = "blast_score_filter_add_gene.csv"

def print_using_time(index_s,start_time):
    # if index_s % 100000 == 0:
    end = time.time() 
    use = end - start_time
    total_line = (index_s + 1) * 100000
    print(f"processed total total_line lines in score file and used time is {use} seconds")
        
def finding_one_score_line_gene(index_s, row_s, sorted_score_df, gtf_df, gtf_index_previous):
     for index_g, row_g in gtf_df[gtf_index_previous:].iterrows():
        # print(index_g,gtf_index_previous)
        # if row_g["source"] == "RefSeq" and row_g["feature"] == "gene":
        start_pos = int(row_s["sstart"]) if int(row_s["sstart"]) <= int(row_s["send"]) else int(row_s["send"])
        end_pos =  int(row_s["send"]) if int(row_s["sstart"]) <= int(row_s["send"]) else int(row_s["sstart"])

        #if score start&end axies in the gene area,then take the gene information to this score line 
        if int(row_g["start"]) <= start_pos and end_pos <= int(row_g["end"]):
            #sorted_score_df.at[index_s,8] = str(row_g[8])
            sorted_score_df.at[index_s, "gene_information"]= row_g["attribute"].replace(',', ';')
            return index_g
        # if score area on the left, gene inoformation area on the right, 
            # means before not in the area,cannot find the gene name, 
            # should return the none found result and also return the inital place for this searching 
        if  start_pos < int(row_g["start"]):
            return gtf_index_previous
    # return gtf_index_previous


def add_gene(input_gtf, input_score, temp_csv, output_csv):
    gtf_df = pd.read_csv(input_gtf, sep="\t", header=None)
    gtf_df.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf_df = gtf_df[(gtf_df["source"] == "RefSeq") & (gtf_df["feature"] == "gene")]
    # Copy the DataFrame before sorting
    gtf_before_sort = gtf_df.copy()
    gtf_df.sort_values(by="start", inplace=True)
    # Check if the DataFrame before and after sorting is the same
    is_same = gtf_df.equals(gtf_before_sort)
    # Print the result
    print("DataFrame the same before and after sorting is", is_same)
    gtf_df.reset_index(drop=True, inplace=True)
    gtf_df["attribute"] = gtf_df["attribute"].str.replace(',', ';')
    print(gtf_df.head())

    temp = pd.read_csv(input_score, sep=",")
    temp.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "qcovhsp", "pident"]
    temp.sort_values(by="sstart", inplace=True)
    temp.reset_index(drop=True, inplace=True)
    temp.to_csv(temp_csv, sep=',', header=True, index=False)
    start_time = time.time() 
    with open(output_csv, 'w') as out_file:
        print_index = 0
        gtf_index = 0
        gtf_index_previous = 0
        for chunk in pd.read_csv(temp_csv, sep=",", chunksize=100000):
            print_using_time(print_index, start_time)
            print_index += 1
            chunk["gene_information"] = np.nan
            chunk["gene_information"] = chunk["gene_information"].astype(str)
            # chunk.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "qcovhsp", "pident"]
            # chunk.sort_values(by="sstart", inplace=True)
            # chunk.reset_index(drop=True, inplace=True)
            # ##add gene column to the sorted_score_df dataframe
            # chunk["gene_information"] = np.nan
            # chunk["gene_information"] = chunk["gene_information"].astype(str)
            # Process each chunk and add gene information

            for index_s, row_s in chunk.iterrows():
                # Process each row in the score file chunk
                gtf_index = finding_one_score_line_gene(index_s, row_s, chunk, gtf_df, gtf_index_previous)
                gtf_index_previous = gtf_index  


            chunk.to_csv(out_file, mode='a', header=None, index=False, quoting=csv.QUOTE_NONE, escapechar='\\')  # Append mode
            
    if os.path.exists(temp_csv):
        os.remove(temp_csv)
        

