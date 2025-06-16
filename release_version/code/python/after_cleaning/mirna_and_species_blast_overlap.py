from io import StringIO
import os
from os import listdir, path, makedirs
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import sys

def compare_lists(list1, list2):
    count = 0
    cnt = 0
    list1 = list1.values.tolist()
    list2 = list2.values.tolist()
    #print(list1,list2,list1[0],list2[0])
    for i in range(0, len(list1)):
        for j in range(cnt, len(list2)):
            if list1[i] < list2[j]:
                break
            if list1[i] == list2[j]:
                cnt = j
                count += 1
    return count

def find_overlap(species_mapping_file, folder_path_microRNA, overlap_file):
    # f = path.join(folder_path_overlap, 'mirna_and_topspecies_overlap_analysis.csv')
    # print(folder_path_microRNA)
    ##microRNA statistics
    files_m = os.listdir(folder_path_microRNA)
    # files_s = os.listdir(folder_path_top_species)
    txt_files_m =  [file for file in files_m if file.endswith('_rna.txt')]
    txt_files_s =  [species_mapping_file]
    txt_files_combined = txt_files_m + txt_files_s
    # print(txt_files_combined)
    total_lines_m = []
    total_lines_s = []
    total_final = []
    file_path = ""
    for file in txt_files_combined:
        #print("one")
        if "rna" in file:
            file_path = os.path.join(folder_path_microRNA, file)
            df = pd.read_csv(file_path, sep='\t', header=None)
        else:
            file_path  = species_mapping_file
            df = pd.read_csv(file_path)
        # df = []
        # # try:
        # df = pd.read_csv(file_path, sep='\t', header= None)
        # except pd.errors.ParserError:
        #     df = pd.read_csv(file_path, sep=',', on_bad_lines='skip', header= None)
            
        df_qseq = df.iloc[:, [0]]
        df_qseq.columns = ['qseqid']

        #print(df_qseq.drop_duplicates().apply(sorted))
        if "rna" in file:
            total_lines_m.append({'unique_qseqied_column':df_qseq.drop_duplicates().apply(sorted), 'unique_qseqid_number_microRNA':df_qseq['qseqid'].nunique(), 'file_name':path.join(folder_path_microRNA, file)})
        else:
            total_lines_s.append({'unique_qseqied_column':df_qseq.drop_duplicates().apply(sorted), 'unique_qseqid_number_species':df_qseq['qseqid'].nunique(), 'file_name':species_mapping_file})
        #print(total_lines_m,total_lines_s)
    for i in range(len(total_lines_m)):
        for j in range(len(total_lines_s)):
            # s_filename = ""
            # if "blast_score" not in total_lines_s[j]['file_name']:
            #     s_filename = total_lines_s[j]['file_name'].replace("blastn_mapping_filter_", "")
            # else:
            s_filename = total_lines_s[j]['file_name']
            m_filename = total_lines_m[i]['file_name']
            # print(s_filename)
            # print(m_filename)
            # print(total_lines_m[i]['unique_qseqied_column'])
            # print(total_lines_s[j]['unique_qseqied_column'])
            # if m_filename == s_filename:
            overlap_count = compare_lists(total_lines_m[i]['unique_qseqied_column'], total_lines_s[j]['unique_qseqied_column'])
            total_final.append({'ovelap-count':overlap_count, 'left-microRNA-count':total_lines_m[i]['unique_qseqid_number_microRNA'] - overlap_count,'right-species-count':total_lines_s[j]['unique_qseqid_number_species'] - overlap_count, 'file_name':"species_mapping_file's path:"+s_filename + " and mirna_mapping_file's path:" + m_filename})

    # print(total_final)
    # Write DataFrame to an Excel file
    total_final = pd.DataFrame(total_final)
    total_final.to_csv(overlap_file, index=False)

def parse_arguments():
    parser = argparse.ArgumentParser(description="overlap the blast results between species mapping and mirna (hairpin) mapping")
    parser.add_argument('-species_mapping_file', required=True,
                        type=str, help='species mappin file (after filtering score and reads length)', default="none")
    parser.add_argument('-hairpin_folder', required=True,
                        type=str, help='hairpin folder', default="none")
    parser.add_argument('-overlap_file', required=True,
                        type=str, help='overlap file path', default="none")                       
    return parser.parse_args()

def main():
    args = parse_arguments()
    if args.species_mapping_file and args.hairpin_folder and args.overlap_file:
        # time_start_s = time.time()
        find_overlap(args.species_mapping_file, args.hairpin_folder, args.overlap_file)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()
    print("finished overlap.")
