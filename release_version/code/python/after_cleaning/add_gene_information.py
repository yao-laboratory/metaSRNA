import re
from io import StringIO
import os
import pandas as pd
import numpy as np
import time
import csv

# import scipy
# from skbio import TreeNode, read

# gcf_file = "gene.gtf"
# score_file = "score_filter.txt"
# file_path = "blast_score_filter_add_gene.csv"

def print_using_time(index_s,start_time):
    if index_s % 100000 == 0:
        end = time.time() 
        use = end-start_time
        print(f"processed {index_s + 1} lines in score file and used time is {use} seconds")
        
def finding_one_score_line_gene(index_s, row_s, sorted_score_df, gcf_df, gcf_index_previous):
    for index_g, row_g in gcf_df[gcf_index_previous:].iterrows():
       # print(index_g,gcf_index_previous)
        if row_g[1] == "RefSeq" and row_g[2] == "gene":
            #if score start&end axies in the gene area,then take the gene information to this score line 
            if int(row_g[3]) <= int(row_s[2]) and int(row_s[3]) <= int(row_g[4]):
                #print(index_g, gcf_index_previous, row_g[3],row_g[4],row_g[8],row_s[2],row_s[3])
                #print(str(row_g[8]))
                #sorted_score_df.at[index_s,8] = str(row_g[8])
                sorted_score_df.at[index_s, 8] = row_g[8]
                #print(sorted_score_df[8])
                return index_g
            # if score area on the left, gene inoformation area on the right, 
              # means before not in the area,cannot find the gene name, 
              # should return the none found result and also return the inital place for this searching 
            if  int(row_s[2]) < int(row_g[3]):
                return gcf_index_previous
    # return gcf_index_previous

def add_gene(input_gcf, input_score, output_csv):
    # gcf_file = "gcf_ecoli.gtf"
    gcf_df = pd.read_csv(input_gcf, sep="\t",header = None)
    # print(gcf_df)
    
    score_df = pd.read_csv(input_score, sep=",",header = None)
    #print(score_df)
    #print(score_df)
    #sortecd with start axis
    # sorted_score_df = score_df[:10000].sort_values(by=2)
    sorted_score_df = score_df.sort_values(by=2)
    #sorted_score_df = score_df.sort_values(by=2, ascending=False).sort_index(level=0, ascending=[True])
    #print(sorted_score_df)
    sorted_score_df.reset_index(drop=True, inplace=True)
    #print(sorted_score_df)
    
    ##add gene column to the sorted_score_df dataframe
    sorted_score_df[8] = np.nan
    sorted_score_df[8] = sorted_score_df[8].astype(str)
    #print(sorted_score_df)
    #set initial timer
    start_time = time.time() 
    #set each time's for loop start position
    gcf_index = 0
    gcf_index_previous = 0
    
    for index_s, row_s in sorted_score_df.iterrows():
        #print time for each 10000 rows processed
        print_using_time(index_s, start_time)
        #finish one score line's gene finding
        gcf_index = finding_one_score_line_gene(index_s, row_s, sorted_score_df, gcf_df, gcf_index_previous)
        gcf_index_previous = gcf_index
    
    #print the final dataframe
    #print(sorted_score_df)
    
    
    #store the final dataframe to txt
    sorted_score_df.to_csv(output_csv, sep=',', header=None, index = None, quoting=csv.QUOTE_NONE, escapechar='\\')

# if __name__ == "__main__":
#     main()
#     print("finished adding.")