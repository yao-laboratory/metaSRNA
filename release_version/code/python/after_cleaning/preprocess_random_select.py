
import random
from Bio import SeqIO
import argparse
def preprocess_random_select(input_fasta, output_fasta, proportion_divisor):

    # First step, count total sequences
    total_sequences = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))
    num_to_select = total_sequences // proportion_divisor  # Calculate 1/4 of total sequences

    # Second step, randomly select sequences
    selected_indices = set(random.sample(range(total_sequences), num_to_select))
    with open(output_fasta, "w") as output_handle:
        for i, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
            if i in selected_indices:
                SeqIO.write(record, output_handle, "fasta")

    print(f"Randomly selected {num_to_select} sequences saved to {output_fasta}.")


def parse_arguments():
    parser = argparse.ArgumentParser(description="random selecte the proportion of clean fasta file as output file")
    parser.add_argument('-input_file', required=True,
                        type=str, help='input file path', default="none")
    parser.add_argument('-output_file', required=True,
                        type=str, help='output file path', default="none")
    parser.add_argument('-select_proportion_divisor', required=True,
                        type=int, help='random select denominator', default="none")

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_file and args.output_file and args.select_proportion_divisor:
        # time_start_s = time.time()
        preprocess_random_select(args.input_file, args.output_file, args.select_proportion_divisor)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()