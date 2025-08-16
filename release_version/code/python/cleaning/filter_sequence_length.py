import argparse
import time

def process_length(input_file, output_folder, minimum_length, maximum_length):
    output1_1 = output_folder + "/final_seq_" + str(minimum_length) + (".fastq" if maximum_length == -1 else "_to_" + str(maximum_length) + ".fastq")
    if maximum_length == -1:
        maximum_length = float('inf')
    with open(input_file, 'r') as file, open(output1_1, 'w') as output_file1_1:
        # Read the first line
        #id1 = 1
        #id2 = 1
        while True:
            lines = []
            # line_unmap = []
            for i in range(4):
                line = file.readline()
                if (not line):
                    break
                # line_unmap.append(line)
                lines.append(line.strip())
                # lines.append(line)
            if (not lines or len(lines) <= 3):
                output_file1_1.close()
                # output_file1_2.close()
                # output_file1_3.close()
                # output_file2_2.close()
                break
            if len(lines) == 4 and "+" in lines[2]:
                # new_string = lines[1].replace("\n", "")
                if len(lines[1]) >= minimum_length and len(lines[1]) <= maximum_length:
                    #line1_w = '@' + str(id1) + '\n'
                    output_file1_1.write(lines[0]+"\n")
                    #id1 += 1
                    output_file1_1.write(lines[1]+"\n")
                    output_file1_1.write(lines[2]+"\n")
                    output_file1_1.write(lines[3]+"\n")

def remove_duplicates(input_file, output_file):
    with open(input_file, 'r') as file, open(output_file, 'w') as output_file1_1:
        sequences = {}
        current_seq = None
        line = file.readline().strip()
        # print(current_value)
        while line and line.startswith('>'):
            # print("a")
            current_key = file.readline().strip()
            if current_key not in sequences:
                sequences[current_key] = line
            line = file.readline().strip()

        # with open(output_file, 'w') as file:
        for key, value in sequences.items():
            output_file1_1.write(f"{value}\n{key}\n")



def main():
    parser = argparse.ArgumentParser(prog='Cleanfastaq')

    subparsers = parser.add_subparsers(
        dest='subcommand', help='Sub Commands Help')
    # add sub command
    parser_c1 = subparsers.add_parser("process_length",
                                      help='')
    parser_c1.add_argument('-input', required=True,
                           type=str, help='input file: fastq format', default="none")
    parser_c1.add_argument('-output_folder', required=True,
                           type=str, help='output file main name', default="none")
    parser_c1.add_argument('-filter_min_length', required=True,
                        type=int, help='minimum length requirement for reads file', default="none")
    parser_c1.add_argument('-filter_max_length', required=False,
                    type=int, help='maximum length requirement for reads file', default=-1)

    # add sub command
    parser_c2 = subparsers.add_parser("remove_duplicates",
                                      help='')
    parser_c2.add_argument('-input', required=True,
                           type=str, help='input file: fastq format', default="none")
    parser_c2.add_argument('-output_file', required=True,
                           type=str, help='output file main name', default="none")                 

    args = parser.parse_args()
    if args.subcommand == 'process_length':
        input_file = args.input
        output_folder = args.output_folder
        min_length = args.filter_min_length
        max_length = args.filter_max_length
        time_start_s = time.time()
        process_length(input_file, output_folder, int(min_length), int(max_length))
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')

    elif args.subcommand == 'remove_duplicates':
        input_file = args.input
        output_file = args.output_file
        time_start_s = time.time()
        remove_duplicates(input_file, output_file)
        time_end_s = time.time()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')
    else:
        print("Wrong input. Check parameters")


if __name__ == "__main__":
    main()