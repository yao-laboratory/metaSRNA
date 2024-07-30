import re
from io import StringIO
import itertools
import sys
import subprocess


def write_one_record(lines, result, output_files, id_map, x_num):
    #print("store_result", result)
    # print("x_num",x_num, len(temp_result))
    for i in range(x_num):
        if i == 0:
            # print( result[0],output_files[0])
            line1_w = '@' + str(id_map) + '\n'
            output_files[0].write(line1_w)
            line2_w = result[0]+'\n'
            output_files[0].write(line2_w)
            line3_w = '+'+'\n'
            output_files[0].write(line3_w)
            filter_q_line = re.sub("^£+", "", lines[3])
            # print("fql", filter_q_line)
            line4_w = filter_q_line[:len(result[0])]+'\n'
            output_files[0].write(line4_w)
            # print(line1_w,line2_w,line3_w,line4_w)
        else:
            # print("i",i)
            line1_w = '@' + str(id_map) + '\n'
            output_files[i].write(line1_w)
            line_special = ''
            if i + 1 > len(result): 
                line_special = '\n'
            else:
                line_special = result[i] + '\n'
            output_files[i].write(line_special)
            

def compare_string(remaining_target_seq,required_target_seq,direction):
    tolerance = 0
    # min_len = min(len(remaining_target_sequences), len(required_target_sequences))
    if len(required_target_seq) <= len(remaining_target_seq):
        for i in range(len(required_target_seq)):
            remaining_target_start_pos = remaining_target_seq[i] if direction == "head" else remaining_target_seq[i-len(required_target_seq)]
            if required_target_seq[i] != remaining_target_start_pos:
                tolerance += 1
    else:
        return -1
    return tolerance


def compare_targets_torlerance(result, pattern, threshold, required_target_seq, information_list):
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num] = information_list
    threshold_minimum = sys.maxsize
    single_final_result = ''
    if result:
        # compare only one combination of regular expression, but there are maybe more than one suitable result(single_result)
        for single_result in result:
            # print("single_result",single_result)
            # print("required_target_seq",required_target_seq)
            # print("match_seq",match_seq)
            threshold_temp = 0
            flag = True
            final_result = list(single_result)
            for i in range(len(required_target_seq)):
                remaining_target_seq = single_result[i+1] if required_target_seq[i][1] == "right" else single_result[i]
                #produce the final string(cut the remaining target sequences) at the same time
                num = len(required_target_seq[i][0])
                if required_target_seq[i][1] == "right":
                    # print("single_result",single_result[i+1][num:])
                    # print("final_result[i]",final_result[i+1])
                    final_result[i+1] =  single_result[i+1][num:]
                else:
                    final_result[i] =  single_result[i][:-num]
                
                direction = "head" if required_target_seq[i][1] == "right" else "tail"
                # there are no possible to have match sequences, cause too short
                if(compare_string(remaining_target_seq,required_target_seq[i][0],direction) == -1):
                    flag = False
                    break
                threshold_temp += compare_string(remaining_target_seq,required_target_seq[i][0],direction)
                #print("threshold temp", threshold_temp)
                if threshold_temp > threshold:
                    break
            if threshold_temp <= threshold and threshold_temp < threshold_minimum and flag == True:
                threshold_minimum = threshold_temp
                single_final_result = final_result
    # print("minimum",threshold_minimum)
    return threshold_minimum, single_final_result


def find_pattern_combination(information_list, pattern):
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num] = information_list
    # print(information_list)
    # find top 5  or last 5 str in one seq string
    options = ['left_characters', 'right_characters']
    # Generate all combinations of choices for each string
    ##this means try to know find the whole pattern or partial pattern
    match_seq_number = seq_num if pattern == "whole_pattern" else (seq_num - 1)
    
    choices = list(itertools.product(options, repeat=match_seq_number)) 
    # print(choices)
    total_combination = []
    target_sequences = []
    for choice in choices:
        concatenated_string = ''
        target_single_sequence = []
        for i in range(match_seq_number):
            
            # calculate how many character being the target for each seq
            cnt = 5 if match_seq_count[i] > 5 else 1
            
            #add first gap requirement before each target sequence
            if i == 0:
                number_part0_x = match_x_count[i] if choice[i] == 'left_characters' else (match_x_count[i] + len(match_seq[i]) - cnt)
                number_part0_star = 1 if choice[i] == 'left_characters' else (1 + len(match_seq[i]) - cnt)
                concatenated_string += r'(\w{' + str(number_part0_star) + r',})' if match_x[i] == '*' else r'(\w{' + str(number_part0_x) + r',})'
           
            #then add regular expression like AACTGw{26}   
            #calculate the gap between the target 5 characters
            number_part1_x = 0
            number_part1_star = 0
            number_1_star =  (len(match_seq[i]) - cnt) if choice[i] == 'left_characters' else 0
            number_1_x = (match_x_count[i+1] + len(match_seq[i]) - cnt) if choice[i] == 'left_characters' else (match_x_count[i+1])
            #judge match_seq already the last one or not, in case beyond the boundary
            if (i + 1) < match_seq_number:
                cnt_next =  5 if match_seq_count[i+1] > 5 else 1
                number_part1_x = number_1_x if choice[i+1] == 'left_characters' else (number_1_x + len(match_seq[i+1]) - cnt_next)
                number_part1_star = number_1_star if choice[i+1] == 'left_characters' else (number_1_star + len(match_seq[i+1]) - cnt_next)
            else:
                number_part1_x = number_1_x
                number_part1_star = number_1_star
            pattern_exactly_count =   r'(\w{' + str(number_part1_star) + r',})'  if match_x[i+1] =='*' else r'(\w{' + str(number_part1_x) + r'})'
            concatenated_string +=  (match_seq[i][:cnt] + pattern_exactly_count) if choice[i] == 'left_characters' else (match_seq[i][-cnt:] + pattern_exactly_count)
            
            ##add left target sequnces to the list
            seq = match_seq[i][cnt:] if choice[i] == 'left_characters' else match_seq[i][:-cnt] 
            pos = "right" if choice[i] == 'left_characters' else "left"
            target_single_sequence.append([seq,pos])
            
        target_sequences.append(target_single_sequence)
        total_combination.append(concatenated_string)
    return target_sequences, total_combination


def match_sequnce_with_torlerance(target_sequences, total_combination, filter_line, information_list, output_files, threshold, id_map):
    threshold_minimum = sys.maxsize 
    final_result = ''
    #make target sequences' rest parts to compare
    for i in range(len(total_combination)):
        pattern = re.compile(total_combination[i])
        result = pattern.findall(filter_line)
        threshold_minimum_temp, final_result_temp = compare_targets_torlerance(result, pattern, int(threshold), target_sequences[i], information_list)
        #if have results, 
        if final_result_temp:
            if final_result_temp[0] is not None and threshold_minimum_temp < threshold_minimum:
                threshold_minimum = threshold_minimum_temp
                final_result = final_result_temp
    return threshold_minimum, final_result


def save_one_record(lines, final_result, threshold_minimum, output_files, id_map, information_list):
    if len(final_result) == 0 or threshold_minimum == sys.maxsize:
        return False
    #save the final strings results to different files
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num] = information_list
    #print("x_num", x_num)
    if match_x[0] != "*":
        #print(final_result[0])
        final_result[0] = final_result[0][:match_x_count[0]]
    # print("final_result", final_result)
    write_one_record(lines, final_result, output_files, id_map, len(match_x))        
    return True 


def find_all_pattern(lines, filter_line, information_list, output_files, threshold, id_map):
    target_sequences, total_combination = find_pattern_combination(information_list, "whole_pattern")
    #print("1,1",target_sequences, total_combination)
    threshold_minimum,final_result = match_sequnce_with_torlerance(target_sequences, total_combination, filter_line, information_list, output_files, threshold, id_map)
    #print("!!!",threshold_minimum,final_result)
    if (save_one_record(lines, final_result, threshold_minimum, output_files, id_map, information_list)):
        return True
    return False
    
        

def find_patial_pattern (lines, filter_line, information_list, output_files, threshold, id_map, tail_torlerance):
    # print("lines",filter_line)
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num] = information_list
    # tail_torlerance = 4
    length = find_unfinished_end_seq(filter_line, match_seq[-1], tail_torlerance)
    #if the final partial target sequence is matched
    if length != -1 and length < len(filter_line):
        #print(length)
        target_sequences, total_combination = find_pattern_combination(information_list,"partial_pattern")
        # print("2,2",target_sequences, total_combination)
        threshold_minimum,final_result = match_sequnce_with_torlerance(target_sequences, total_combination, filter_line[:-length], information_list, output_files, threshold, id_map)
        #print("!!!2",threshold_minimum,final_result)
        if (save_one_record(lines, final_result, threshold_minimum, output_files, id_map, information_list)):
            return True
    return False
    

def find_unfinished_end_seq(actual_seq, required_seq, tail_torlerance):
    length = 0
    # print(sequence, pattern)
    #print("required_seq",required_seq)
    for i in range(len(required_seq), 0, -1):
        # print(pattern[:i],i)
        if (actual_seq.endswith(required_seq[:i])):
            length = i
            break
    
    # length = len(actual_seq) - num
    # print(length,num)
    if length < tail_torlerance:
        length = -1
        
    #print("l",length )
    return length
 
    
    
def find_id(output_file_step1):
    with open(output_file_step1, 'r') as file:
        last_four_lines = file.readlines()[-4:]
        for line in last_four_lines:
            if line.startswith('@'):
                line = line.strip()
                match = re.search(r'\d+', line)
                if match:
                    return int(match.group()) + 1
    return 1



def store_unmatched_file(lines, unmatched_file_step2):
    for i in range(len(lines)):
        unmatched_file_step2.write(lines[i] + '\n')
        

def clean_step2(information_list, input_file, output_file_step1, output_files_step2, unmatched_file_step2, threshold, tail_tolerance):
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num]  = information_list
    with open(input_file, 'r') as file:
        # Read the first line
        id_map = find_id(output_file_step1)
        # id_unmap = 1
        while True:
            lines = []
            line = file.readline()
            #if the file end in this line, then break the loop
            if not line:
                break
            
            if "@" in line:
                lines.append(line.strip())
                for i in range(3):
                    line = file.readline()
                    if (not line):
                        break
                    lines.append(line.strip())
                if (not lines or len(lines) < 4):
                    break

                if len(lines) == 4 and "@" in lines[0] and "+" in lines[2]:
                    filter_line = re.sub("^N+|\s", "", lines[1])

                    if (find_all_pattern(lines, filter_line, information_list, output_files_step2, threshold, id_map)):
                        id_map += 1
                    else:
                        if (find_patial_pattern (lines, filter_line, information_list, output_files_step2, threshold, id_map, tail_tolerance)):
                            id_map += 1
                        else:
                            store_unmatched_file(lines, unmatched_file_step2)
    file.close()




