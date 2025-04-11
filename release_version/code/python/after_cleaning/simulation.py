import matplotlib.pyplot as plt
from PIL import Image
from brokenaxes import BrokenAxes
from intervaltree import Interval, IntervalTree
import numpy as np
import os
import random
import matplotlib.cm as cm
import argparse

def merge_intervals(intervals):
    """Merge overlapping x-axis intervals using IntervalTree."""
    tree = IntervalTree(Interval(start, end) for start, end in intervals)
    #merge overlapping intervals
    tree.merge_overlaps() 
    #return list (not tuple) 
    return sorted((iv.begin, iv.end) for iv in tree)

# function to plot each dataset
def plot_dataset(dataset, dataset_idx, output_folder):
    num_blocks = len(dataset)
    cmap = cm.get_cmap('tab20', num_blocks)
    fig, ax = plt.subplots(figsize=(12, 4))
    for block_index, block in enumerate(dataset):
        x_min = min(int(item[0]) for item in block)
        x_max = max(int(item[1]) for item in block)
        #print(f"Dataset {dataset_idx}: xlim = ({x_min}, {x_max})")
        parameter = 1
        # Step 1: Create unique x-limits based on dataset
        important_ranges = []
        important_ranges = [(start, end) for start, end, _ in block]
        merged_xlims = merge_intervals(important_ranges)
        print(merged_xlims)
        # Step 3: Plot horizontal lines
        distance = 0
        color = cmap(block_index)
        for start, end, line_count in block:
            print(start,end,line_count)
            if merged_xlims and int(start) == merged_xlims[0][0]:
                # logging.info("start")
                # logging.info(start)
                distance = 0
                merged_xlims.pop(0)
            y_positions = [distance + parameter * (i + 1) for i in range(int(line_count))]
            print(y_positions)
            # random_color = np.random.rand(3,)
            for y in y_positions:
                plt.hlines(y, start, end, colors=color, linestyles='solid', linewidth=1)
            distance = y_positions[-1]
            

    ax.set_xlabel("Genome Position")
    ax.set_ylabel("Coverage")    
    # plt.show()

    #save the plot to a buffer as a grayscale image
    buffer_path = f"{output_folder}/simulation_1_200_gaps_{dataset_idx}_blocks.png"
    plt.savefig(buffer_path, format='png', dpi=300, bbox_inches='tight')
    plt.close()

def collect_all_blocks_data(simulation_data_folder, files):
    bed_blocks = {}
    block_lines = []
    print(files)
    for file_name in files:
        try:
            file_path = os.path.join(simulation_data_folder, file_name)
            with open(file_path, "r") as file:
                bed_blocks[file_name] = []
                chrom_blocks = []
                line_block = []
                collecting = False
                for line in file:
                    line = line.strip()
                    if line == "":
                        collecting = False
                        if line_block:
                            bed_blocks[file_name].append(line_block)
                            line_block = []
                    elif line == "original bed lines:":
                        collecting = True
                    elif collecting:
                        columns = line.split('\t')
                        line_block.append([columns[0],columns[1],columns[2],columns[-1]])
        except Exception as e:
            print(f"Error opening file '{file_paths.get(file_name, 'UNKNOWN')}': {e}")
    return bed_blocks

def cluster_blocks_by_chrom(bed_blocks, files):
    group_all = {}
    for file_name in files:
        group = {}
        # group[chrom] = []
        group_all[file_name] = {}
        for block in bed_blocks[file_name]:
            chrom = block[0][0]
            if chrom not in group:
                group[chrom] = []
            group[chrom].append(block)
        group_all[file_name] = group
    return group_all

def simulate_data_generation(original_data, blocks_number, start_point):
    random_block_id = random.randint(0, blocks_number - 1)
    original_min_start_point = min(int(line[1]) for line in original_data[random_block_id])
    original_max_end_point = max(int(line[2]) for line in original_data[random_block_id])
    move_units = original_min_start_point - start_point
    max_end_point = original_max_end_point - move_units
    block_data = []
    for line in original_data[random_block_id]:
        line_start = int(line[1])
        line_end = int(line[2])
        line_count = int(line[3])
        # print(line[1],line[2],line[3],min_start_point)
        count = 50 if line_count > 50 else line_count
        # print(line_start - min_start_point, line_end - min_start_point, count)
        block_data.append([line_start - move_units, line_end - move_units, count])
    
    return block_data, max_end_point

def generate_all_data(bed_blocks_group_by_chrom, files, blocks_number, smallest_gap, biggest_gap):
    number = 0
    dataset = []
    random_file = random.choice(files)
    random_chrom = random.choice(list(bed_blocks_group_by_chrom[random_file].keys()))
    start_point = 1 
    while number < blocks_number:
        random_file = random.choice(files)
        actual_blocks_number = len(bed_blocks_group_by_chrom[random_file][random_chrom])
        block_data, end_point = simulate_data_generation(bed_blocks_group_by_chrom[random_file][random_chrom], actual_blocks_number, start_point)
        dataset.append(block_data)
        start_point = end_point + random.randint(smallest_gap, biggest_gap)
        number += 1
    return dataset


def simulate_blocks(input_files, number_of_blocks, block_gap_min_threshold, block_gap_max_threshold, output_folder):
    # simulation_data_folder = "simulation_data"
    # files = ["fully_symmetric_blocks.txt", "fully_overlapping_blocks.txt", 
    #             "partially_overlapping_blocks.txt", "singleton_blocks.txt", "other_situation_blocks.txt"]
    #prepare all the data
    input_files_list = input_files.split()
    if not input_files_list:
        raise ValueError(f"input files' paths: '{input_files}' is wrong")
    files = []
    for file_path in input_files_list:
        files.append(os.path.basename(file_path))
    simulation_data_folder = os.path.dirname(input_files_list[0])

    bed_blocks = collect_all_blocks_data(simulation_data_folder, files)
    # print(bed_blocks[files[0]])
    bed_blocks_group_by_chrom = cluster_blocks_by_chrom(bed_blocks, files)
    # print(bed_blocks_group_by_chrom)
    # print(bed_blocks_group_by_chrom[files[0]])
    #start simulating data
    simulation_dataset = generate_all_data(bed_blocks_group_by_chrom, files, int(number_of_blocks), int(block_gap_min_threshold), int(block_gap_max_threshold))
    print(simulation_dataset)
    plot_dataset(simulation_dataset, 15, output_folder)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Simulate the real genome mapping data and segment it into blocks by color using input files from various block classification categories. The output is a PNG image.")
    parser.add_argument('-input_files', required=False,
                        type=str, help="Input files' paths sperate by space ")
    parser.add_argument('-number_of_blocks', required=False,
                        type=str, help='Decide how many blocks to simulate.')
    parser.add_argument('-block_gap_min_threshold', required=False,
                        type=str, help='Min threshold about the gaps between the blocks.')
    parser.add_argument('-block_gap_max_threshold', required=False,
                        type=str, help='Max threshold about the gaps between the blocks.')
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output resutls folder path')

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_files and args.number_of_blocks and args.block_gap_min_threshold and args.block_gap_max_threshold and args.output_folder:
        # time_start_s = time.time()
        simulate_blocks(args.input_files, args.number_of_blocks, args.block_gap_min_threshold, args.block_gap_max_threshold, args.output_folder)
        # time_end_s = time.time()
        # time_c = time_end_s - time_start_s
        # print('time cost', time_c, 's')


if __name__ == "__main__":
    main()