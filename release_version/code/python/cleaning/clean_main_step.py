import argparse
import time
import re
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from clean_fastq_step1 import clean_step1
from clean_fastq_step2 import clean_step2

# def find_start_pos(string):
#     new_string = [] 
#     for i in range(len(string)):
#         if i == 0:
#             new_string.append(string[i])
#         if i + 1 >= len(string):
#             break
#         elif string[i+1] == string[i] + 1:
#             i += 1
#             continue
#         else:
#             new_string.append(string[i + 1])
#     return new_string
            

def analyze_format(fa_format):
    # pattern = r'[xX]'
    # p_x = [(match.start()) for match in re.finditer(pattern, fa_format)]
    # pattern = r'[*]'
    # p_star = [(match.start()) for match in re.finditer(pattern, fa_format)]
    # pattern = r'[ATCGatcg]'
    # p_seq = [(match.start()) for match in re.finditer(pattern, fa_format)]
    # print("!!!!!!!!")
    # print("p_x", p_x)
    # x_count = find_start_pos(p_x)
    # star_count = find_start_pos(p_star)
    # seq_count = find_start_pos(p_seq)
    # Convert pattern's small letter to big letter
    fa_format = re.sub(r'x', 'X', fa_format)
    fa_format = re.sub(r'a', 'A', fa_format)
    fa_format = re.sub(r't', 'T', fa_format)
    fa_format = re.sub(r'c', 'C', fa_format)
    fa_format = re.sub(r'g', 'G', fa_format)

    match_x = re.findall(r'[X*]+', fa_format)
    # match_star = re.findall(r'[*]+', fa_format)
    match_seq = re.findall(r'[ATCG]+', fa_format)
    match_x_count = [len(s) for s in match_x]
    # star_count = [len(s) for s in match_star]
    match_seq_count =[len(s) for s in match_seq]
    # file_x_num = 0
    # file_start_num = 0
    # if match_x == None:
    #     file_x_num = 0
    # else:
    #     file_x_num = len(match_x)
    # if match_star == None:
    #     file_star_num = 0
    # else:
    #     file_star_num = len(match_star)
    # file_num = len(match_x) + len(match_star)
    # if match_seq == None:
    #     file_seq_num = 0
    # else:
    #     file_seq_num = len(match_seq)
    x_num = len(match_x_count)
    seq_num = len(match_seq_count)
    
    # print(p_x)
    # print(p_star)
    # print(p_seq)
    # print(start_x)
    # print(start_star)
    # print(start_seq)
    # print(match_x)
    # print(match_star)
    # print(match_seq)
    print([match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num])
    return [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num]

def clean_start(input_file,fa_format):
    information_list = analyze_format(fa_format)
    file_num = information_list[-2]
    num_files = int(information_list[-2]) if len(information_list) == 6 and str(information_list[-2]).isdigit() else 0
    seq_num = information_list[-1]
    if (file_num == None) or (seq_num == None) or (num_files == 0):
        print("input sequence formatting wrong")
        return False
    # f = 'MEE-OMVs-1_S4_L001_R1_001.fastq'
    # output1_1 = "output_cleaned_seq.fastq"
    # output1_2 = "output_umi_seq.fastq"
    # output1_3 = "output_unmatched_seq.fastq"
    # output2_1 = "output_seq_and_umi12.csv"
    # information_list = [start_x,match_x,start_star,match_star,start_seq,match_seq,file_num,seq_num]
    #print(information_list)
    return information_list

def clean_main_step(input_file,output_filename,fa_format, tolerance, tail_tolerance):
    information_list = clean_start(input_file, fa_format)
    if (information_list):
        # output_files = []
        #find create a new folder
        output_files_step1 = [open(output_filename + f"_{i}_step1.fastq",'w') for i in range(1, int(information_list[-2])+1)]
        if any('*' in item and 'X' in item for item in information_list[1]):
            output_files_step1.append(open(output_filename + f"_step1_umi.fastq",'w'))
        unmatched_file_step1 = open(output_filename + f"_step1_unmatched.fastq",'w') 
        # print("2",unmatched_files)
        # sorted(output_files)
        # output_files.append(unmatched_file)
        # print("1",output_files)
        clean_step1(information_list, input_file, output_files_step1, unmatched_file_step1, tail_tolerance)
        #close the files first then open in next step
        for i in range(len(output_files_step1)):
            output_files_step1[i].close()
        unmatched_file_step1.close()
        print("cleaning_step 1 end.")
        if (int(tolerance) > 0):
            output_files_step2 = [open(output_filename + f"_{i}_step2.fastq",'w') for i in range(1, int(information_list[-2])+1)]
            if any('*' in item and 'X' in item for item in information_list[1]):
                output_files_step2.append(open(output_filename + f"_step2_umi.fastq",'w'))
            unmatched_file_step2 = open(output_filename + f"_unmatched.fastq",'w') 
            input_file = output_filename + f"_step1_unmatched.fastq"
            output_step1_firstfile = output_filename + f"_1_step1.fastq"
            clean_step2(information_list, input_file, output_step1_firstfile, output_files_step2, unmatched_file_step2, tolerance, tail_tolerance)
            #close the files    
            for i in range(len(output_files_step2)):
                output_files_step2[i].close()
            unmatched_file_step2.close()
            print("cleaning_step 2 finished.")

def after_clean_step(input_file, output_path):
    with open(input_file, 'r') as file, open(output_path, 'w') as output_file:
        # Read the first line
        id1 = 0
        while True:
            lines = []
            line = file.readline()
            if not line:
                break
            if "@" in line:
                lines.append(line.strip())
                for i in range(3):
                    line = file.readline()
                    #if the file end in this line, then break the loop
                    if not line:
                        break
                    lines.append(line.strip())
                id1 += 1
                if len(lines) == 4 and "@" in lines[0] and "+" in lines[2]:
                    line1_w = '@' + str(id1) + '\n'
                    output_file.write(line1_w)
                    line2_w = lines[1]+'\n'
                    output_file.write(line2_w)
                    line3_w = '+'+'\n'
                    output_file.write(line3_w)
                    line4_w = lines[3]+'\n'
                    output_file.write(line4_w)

def clean_N_seqs(input_file, output_path):
    with open(output_path, "w") as out_file:
        new_index = 1  
        # Iterate through all sequences in the input fastq file
        for record in SeqIO.parse(input_file, "fastq"):
            sequence = str(record.seq).replace('n', 'N')
            ##drop the start and end Ns in the sequences
            trimmed_sequence = sequence.strip("N")
            ##if N in the middle of the sequence, directly drop
            if len(trimmed_sequence) == 0 or "N" in trimmed_sequence:
                continue

            new_record = SeqRecord(
                Seq(trimmed_sequence),             
                id=str(new_index),                         
                description="",
                letter_annotations={"phred_quality": record.letter_annotations["phred_quality"][(len(record.seq) - len(record.seq.lstrip("N"))):len(record.seq.rstrip("N"))]}
            )
            # Write the updated record to the output file
            SeqIO.write(new_record, out_file, "fastq")
            
            # Increment the new index
            new_index += 1

    print(f"Processing complete. clean N in original FASTQ file saved as '{output_path}'.")
