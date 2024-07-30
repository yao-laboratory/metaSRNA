import re
from io import StringIO

#import scipy
#from skbio import TreeNode, read

def write_fastq(lines, temp_result, output_files, id, x_num):
    # print("store_results",temp_result)
    # print("x_num",x_num, len(temp_result))
    for i in range(x_num):
        if i == 0:
            line1_w = '@' + str(id) + '\n'
            output_files[0].write(line1_w)
            line2_w = temp_result[0]+'\n'
            output_files[0].write(line2_w)
            line3_w = '+'+'\n'
            output_files[0].write(line3_w)
            filter_q_line = re.sub("^£+", "", lines[3])
            # print("fql", filter_q_line)
            line4_w = filter_q_line[:len(temp_result[0])]+'\n'
            output_files[0].write(line4_w)
        else:
            # print("i",i)
            line1_w = '@' + str(id) + '\n'
            output_files[i].write(line1_w)
            line_special = ''
            if i + 1 > len(temp_result): 
                line_special = '\n'
            else:
                line_special = temp_result[i] + '\n'
            output_files[i].write(line_special)
            
        
# def write_part1_fastq(temp_result, output_files):
#     line_special = temp_result[1] + '\n'
#     output_files[1].write(line_special)
    
# def write_part2_fastq(temp_result, output_files):
#     line_end = temp_result[2] + '\n'
#     output_files[2].write(line_end)

def meet_formatting_requirement(lines, output_files, result, information_list, id):
    # print("needed",lines, output_files, result, file_num)
    # print("original_line",lines[1])
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num] = information_list
    # print("information_list",[match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num])
    temp_result = []
    # print("result[0]",result[0])
    # print("len(result)",len(result))
    if seq_num == 1:
        for i in range(x_num):
            num = match_x_count[i]
            # print(num,result[i])
            if '*' not in match_x[i] :
                if len(result[i]) < num:
                    return False
                else:
                    temp_result.append(result[i][:num])
            else:
                temp_result.append(result[i])   
    elif seq_num == 2:
        for i in range(x_num):
            num = match_x_count[i]
            if i == 1:
                if '*' not in match_x[i] :
                    if len(result[i]) != num:
                        return False
                    else:
                        temp_result.append(result[i])
                else:
                    temp_result.append(result[i]) 
            else:                
                if '*' not in match_x[i] :
                    if len(result[i]) < num:
                        return False
                    else:
                        temp_result.append(result[i][:num])
                else:
                    temp_result.append(result[i])   
    elif seq_num == 3:
         for i in range(x_num):
            num = match_x_count[i]
            if i == 1 or i == 2:
                if '*' not in match_x[i] :
                    if len(result[i]) != num:
                        return False
                    else:
                        temp_result.append(result[i])
                else:
                    temp_result.append(result[i]) 
            else:
                if '*' not in match_x[i] :
                    if len(result[i]) < num:
                        return False
                    else:
                        temp_result.append(result[i][:num])
                else:
                    temp_result.append(result[i])   
    
    write_fastq(lines, temp_result, output_files, id, len(match_x))
    return True

def find_all_correct_pattern(lines, p_string, filter_line, information_list, output_files, id):
    pattern = re.compile(p_string)
    result = pattern.findall(filter_line)
    # print("p_string",p_string,filter_line, id2)
    #print("1",pattern,result)
    if (result and result[0]):
        #print("2")
        [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num] = information_list
        # print("result",result)
        # print("result_0",result[0])
        if (meet_formatting_requirement(lines, output_files, result[0], information_list, id)):
            return True
    return False
    
def find_pattern(sequence, pattern, min_tolerance):
    num = 0
    # print(sequence, pattern)
    for i in range(len(pattern), 0, -1):
        # print(pattern[:i],i)
        if sequence.endswith(pattern[:i]):
            num = i
            break
    length = len(sequence) - num
    # print(length,num)
    if num + 1 < min_tolerance:
        length = -1
    return length
 
def find_patial_correct_pattern(lines, p_string_part, match_seq, filter_line, information_list,
                                output_files, id, min_tolerance):
    # print("partial start")
    
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num]  = information_list
    final_part_sequence = ''
    final_result = ()
    # > 1 required sequence
    if "w*" in p_string_part:
        temp_list = []
        pattern = re.compile(p_string_part)
        result = pattern.findall(filter_line)
        # print(lines,p_string,p_string_part)
        if result and result[0]: 
            final_part_sequence = result[0][-1]
            final_string_length = find_pattern(final_part_sequence, match_seq[-1], min_tolerance)
            if final_string_length != -1:
                list_final_string = list(result[0])
                list_final_string[-1] = list_final_string[-1][:final_string_length]
                final_result = tuple(list_final_string)
    # only 1 required sequence
    else:
        final_part_sequence = p_string_part
        # print(final_part_sequence,match_seq[-1])
        final_string_length = find_pattern(final_part_sequence, match_seq[-1], min_tolerance)
        # print(final_string_length)
        if final_string_length != -1:
            final_part_sequence = final_part_sequence[:final_string_length]
            final_result = ()
            final_result = final_result + (final_part_sequence,)
    if len(final_result) > 0:
        information_list = [match_x_count, match_x, match_seq_count, match_seq, seq_num, seq_num]
        if (meet_formatting_requirement(lines, output_files, final_result, information_list, id)):
            return True
    return False


def store_unmatched_file(lines, unmatched_file_step1):
    for i in range(len(lines)):
        unmatched_file_step1.write(lines[i] + '\n')
        

def clean_step1(information_list, input_file, output_files_step1, unmatched_file_step1, min_tolerance):
    [match_x_count, match_x, match_seq_count, match_seq, x_num, seq_num]  = information_list
    #print(information_list)
    with open(input_file, 'r') as file:
        # Read the first line
        id1 = 1
        id2 = 1
        while True:
            lines = []
            # line_unmap = []
            line = file.readline()
            #if the file end in this line, then break the loop
            #print(line)
            if not line:
                # print("id1",id1)
                # print("id2",id2)
                break
            if "@" in line:
                #print("id2",id2)
                lines.append(line.strip())
                for i in range(3):
                    line = file.readline()
                    #if the file end in this line, then break the loop
                    if not line:
                        # print("id1",id1)
                        break
                    # line_unmap.append(line)
                    lines.append(line.strip())
                id2 +=1
                    # lines.append(line)
                # if (not lines or len(lines) <= 3):
                #     # output_file1_1.close()
                #     # output_file1_2.close()
                #     # output_file1_3.close()
                #     # output_file2_2.close()
                #     break
                # print("id2",id2,lines,len(lines),lines[0],lines[2])
                if len(lines) == 4 and "@" in lines[0] and "+" in lines[2]:
                    filter_line = re.sub("^N+|\s", "", lines[1])
                    p_string = ""
                    for i in range(len(match_seq)):
                        p_string += r'(\w*)' + match_seq[i]
                    p_string += r'(\w*)'
                    # print(p_string)
                    # print(p_string, id2)
                    if (find_all_correct_pattern(lines, p_string, filter_line, information_list, output_files_step1, 
                                                 id1)):
                        id1 += 1
                    else:
                        if seq_num == 1:
                            p_string_part = filter_line
                        else:
                            p_string_part = ""
                            for i in range(len(match_seq) - 1):
                                p_string_part += r'(\w*)' + match_seq[i]
                            p_string_part += r'(\w*)'
                            # print(p_string_part)

                        if (find_patial_correct_pattern(lines, p_string_part, match_seq, filter_line,     
                                                        information_list, output_files_step1, id1, min_tolerance)):
                            id1 += 1
                        else:
                            store_unmatched_file(lines, unmatched_file_step1)
    file.close()


          


