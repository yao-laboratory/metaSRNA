import os
from os import listdir, path, makedirs
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
# import Levenshtein
from Bio import Phylo, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Clade
from Bio.SeqUtils import seq3
from skbio import TreeNode
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import pandas as pd
import shutil
import sys
# import time
from datetime import datetime
from multiprocessing import Pool, cpu_count
from itertools import combinations
from scipy.sparse import lil_matrix

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from PIL import Image

from sklearn.metrics.pairwise import euclidean_distances
from tensorflow.keras.applications import VGG16
from tensorflow.keras.preprocessing import image
from tensorflow.keras.applications.vgg16 import preprocess_input
from sklearn.metrics import silhouette_score
from intervaltree import Interval, IntervalTree
import pandas as pd
import logging
# from keras.applications.vgg16 import VGG16, preprocess_input
# from keras.preprocessing import image
##set tree depth limit is 5000
sys.setrecursionlimit(5000)
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
# Define the sliding window length
SLIDING_WINDOW_LENGTH = 10000
ABNORMAL_TOLERANCE_FACTOR = 10
NATURAL_GAP_DISTANCE = 300
OVERLAP_BETWEEN_WINDOWS = 120

def process_bed_file_in_ranges(file_path):
    # read the BED file as a DataFrame
    bed_whole_df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'sequence', 'id', 'ids', 'same_seq_count'])
    ##iterate each chrom data
    for chromsome in bed_whole_df['chrom'].unique():
        bed_df = bed_whole_df[bed_whole_df['chrom'] == chromsome]
        min_start_pos = bed_df['start'].min()
        max_end_pos = bed_df['end'].max()

        # the number of rows in the DataFrame
        num_rows = bed_df.shape[0] 
        # the number of sliding windows
        num_sliding_windows = (max_end_pos - min_start_pos) // SLIDING_WINDOW_LENGTH + 1
        # average number of lines per window
        avg_lines_per_window = num_rows // num_sliding_windows + 1
        # maximum tolerance for line numbers
        max_tolerance_lines = avg_lines_per_window * ABNORMAL_TOLERANCE_FACTOR
        logging.info(f"max tolerance lines number are {max_tolerance_lines} ")
        current_start = min_start_pos
        ##deal with each sliding window's data
        while current_start < max_end_pos:
            logging.info("start a new sliding window")
            current_end = current_start + SLIDING_WINDOW_LENGTH
            max_prev_end = 0
            range_df_temp = bed_df[(bed_df['start'] >= current_start) & (bed_df['end'] < current_end)]
            # Keep track of the previous row's end
            prev_end = 0  
            temp_positions = []
            temp_sequences = []
            temp_ids = []
            range_df = pd.DataFrame()
            if range_df_temp.shape[0] > max_tolerance_lines:
                len_df_before = len(range_df_temp)
                logging.info(f"before sampling: {len_df_before}")
                range_df = range_df_temp.sample(n=min(max_tolerance_lines, len(range_df_temp)), random_state=42)
                len_df_after = len(range_df)
                logging.info(f"after sampling: {len_df_after}")
            else:
                range_df = range_df_temp
            range_df = range_df.sort_values(by="start")
            range_df = range_df.reset_index(drop=True)
            return_index = 0
            #this for loop for each sliding window
            for idx, row in range_df.iterrows():
                ### seperate to different block if has natual gap
                if max_prev_end != 0:
                    ###seperate blocks from big natual gap
                    if (int(row['start']) - NATURAL_GAP_DISTANCE) > (int(max_prev_end)):
                        logging.info(f"reblock when meet {NATURAL_GAP_DISTANCE} gap")
                        return_df = range_df.loc[return_index:idx-1]
                        yield temp_ids, temp_sequences, temp_positions, return_df
                        # reset collections
                        temp_sequences = []
                        temp_ids = []
                        temp_positions = []
                        return_index = idx

                temp_positions.append((int(row['start']), int(row['end'])))
                temp_sequences.append(row['sequence'])
                temp_ids.append(row['id'])
                
                max_prev_end = max(max_prev_end, int(row['end']))
                

            # yield the rest batch
            if temp_positions:
                yield temp_ids, temp_sequences, temp_positions, range_df.loc[return_index:]
            
            # Move to the next range
            current_start += SLIDING_WINDOW_LENGTH - OVERLAP_BETWEEN_WINDOWS

# def process_bed_file_in_ranges(file_path, start=0, overlap=120):
#     # Read the BED file as a DataFrame
#     bed_df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'sequence', 'id', 'ids', 'same_seq_count']) 
#     # print("bed_df", bed_df)
#     # Get the maximum end position to define loop termination
#     max_position = bed_df['end'].max()  
#     # Process the BED file in specified ranges
#     current_start = start

#     while current_start < max_position:
#         current_end = current_start + SLIDING_WINDOW_LENGTH
#         # Filter the DataFrame to include only rows in the current range
#         range_df = bed_df[(bed_df['start'] >= current_start) & (bed_df['end'] < current_end)]
#         # print("range_df",range_df)
#         positions = []
#         sequences = []
#         ids = []
#         # same_seq_count = []
#         for _, row in range_df.iterrows():
#             positions.append((int(row['start']), int(row['end'])))
#             sequences.append(row['sequence'])
#             ids.append(row['id'])
#             # same_seq_count.append(row['same_seq_count'])
#         # Return the results for the current range
#         yield ids, sequences, positions, range_df
        
#         # Move to the next range
#         current_start += SLIDING_WINDOW_LENGTH - overlap


# # read bed file
# def read_bed_file(file_path):
#     ids = []
#     sequences = []
#     positions = []
#     same_seq_count = []
#     with open(file_path, 'r') as file:
#         for line in file:
#             if line.strip():
#                 parts = line.strip().split()
#                 if len(parts) >= 4:
#                     start = int(parts[1])
#                     end = int(parts[2])
#                     positions.append((start, end))
#                     sequences.append(parts[3])
#                     ids.append(parts[4])
#                     same_seq_count.append(parts[6])
#     return ids, sequences, positions

def vectorized_overlap_check(positions, output_folder, num_seqs):
    # Convert positions to NumPy arrays for vectorized operations
    positions = np.array(positions)
    # Extract start positions
    starts = positions[:, 0]  
    # Extract end positions
    ends = positions[:, 1]   
    # Compute overlap matrix using broadcasting
    overlap_matrix = (np.maximum.outer(starts, starts) <= np.minimum.outer(ends, ends))

    return overlap_matrix

def vetorized_sequences_gap(positions):
    # convert positions to NumPy arrays
    positions = np.array(positions)
    # positions = np.array(positions)
    starts = positions[:, 0]
    ends = positions[:, 1]  
    # print("starts", starts)
    # print("ends", ends)
    # Compute the gap: start_j - end_i
    gap_matrix = np.subtract.outer(starts, ends)  
    gap_matrix = np.where(gap_matrix > 0, gap_matrix, 0)  
    # add lower triangle values to the upper triangle for symmetry
    gap_matrix = gap_matrix + gap_matrix.T

    return gap_matrix

# function to align and check compensatory match
def align_and_check(params):
    i, j, seq_i_rc, seq_j, threshold, aligner = params
    # if distance_between_sequences >= 100:
    #     return i, j, False
    alignments = aligner.align(seq_i_rc, seq_j)
    if not alignments:
        return i, j, False

    best_alignment = alignments[0]
    aligned_seq1 = best_alignment.aligned[0]
    aligned_seq2 = best_alignment.aligned[1]

    match_count = sum(
        seq_i_rc[a] == seq_j[b]
        for block1, block2 in zip(aligned_seq1, aligned_seq2)
        for a, b in zip(range(block1[0], block1[1]), range(block2[0], block2[1]))
    )
    shorter_length = min(len(seq_i_rc), len(seq_j))
    is_compensatory = (match_count / shorter_length) >= threshold

    return i, j, is_compensatory

# chunk generator to avoid large memory usage
def generate_chunks(pos_distances, num_seqs, sequences, threshold, aligner, chunk_size):
    combinations_generator = combinations(range(num_seqs), 2)
    chunk = []
    for i, j in combinations_generator:
        seq_i_rc = str(Seq(sequences[i]).reverse_complement())
        seq_j = sequences[j]
        sequences_gap = pos_distances[i,j]
        if sequences_gap > 2 and sequences_gap <= 100:
            chunk.append((i, j, seq_i_rc, seq_j, threshold, aligner))
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
        
def vectorized_compensatory_check_parallel(sequences, pos_distances, threshold=0.7):
    num_seqs = len(sequences)

    # initialize sparse compensatory matrix
    compensatory_matrix = lil_matrix((num_seqs, num_seqs), dtype=bool)

    # initialize aligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    # chunked task generation and parallel processing
    chunk_size = 1000  
    with Pool(cpu_count()) as pool:
        for chunk in generate_chunks(pos_distances, num_seqs, sequences, threshold, aligner, chunk_size):
            results = pool.map(align_and_check, chunk)
            for i, j, is_compensatory in results:
                if is_compensatory:
                    compensatory_matrix[i, j] = compensatory_matrix[j, i] = True

    return compensatory_matrix

def get_bed_lines_by_indices(bedfile_partial_df, indices):
    # # read all lines in bed
    # with open(bed_file_path, 'r') as file:
    #     bed_lines = [line.strip() for line in file if line.strip()]
    # print("bedfile_partial_df before:",bedfile_partial_df)
    # error

    if max(indices) >= len(bedfile_partial_df):
        raise IndexError("index out of the bed lines")

    #get the lines according ids index 
    # selected_lines = [bedfile_partial_df[i] for i in indices]
    # selected_lines = bedfile_partial_df.iloc[indices]
    # logging.info("indices:")
    # logging.info(indices)
    # logging.info(selected_lines)
    sorted_indices = sorted(indices)
    # logging.info("sorted_indices:")
    # logging.info(sorted_indices)
    selected_lines = bedfile_partial_df.iloc[sorted_indices]
    # selected_lines.sort_values(by=['chrom', 'start', 'end'], ascending=[True, True, True], inplace=True)
    # logging.info("sorted selected_lines:")
    logging.info(selected_lines)
    # selected_lines = bedfile_partial_df.iloc[indices]
    return selected_lines

def add_qseqid_to_list(qseqid_list,member_ids):
    for id in member_ids:
        if int(id) not in qseqid_list:
            qseqid_list.append(int(id))

def compute_distance_matrix(sequences, positions, output_folder):
    num_seqs = len(sequences)
    if num_seqs == 1:
        logging.info("only one sequence available, return zero matrix.")
        return np.zeros((1, 1)), np.zeros((1, 1), dtype=bool), np.zeros((1, 1), dtype=bool)
    dist_matrix = np.zeros((num_seqs, num_seqs))
    #dist_matrix = np.memmap(path.join(output_folder,"middle_results/dist_matrix.dat"), dtype=np.int32, mode='w+', shape=(num_seqs, num_seqs))
    ###old code start 
    overlap_matrix = vectorized_overlap_check(positions, output_folder, num_seqs)
    ###old code end
    # overlap_matrix = np.memmap('overlap_matrix.dat', dtype=bool, mode='w+', shape=(num_seqs, num_seqs))

    #print("finished calculating distance 1 step",datetime.now() )
    pos_distances_matrix = vetorized_sequences_gap(positions)
    #print("finished calculating distance 2 step",datetime.now() )
    compensatory_matrix = vectorized_compensatory_check_parallel(sequences, pos_distances_matrix)
    #print("finished calculating distance 3 step",datetime.now() )
    compensatory_matrix_dense = compensatory_matrix.toarray()

    # compute positional distances
    # compute middle positions
    mid_positions = np.array([(pos[0] + pos[1]) / 2 for pos in positions])
    #print("finished calculating distance 3 step",datetime.now() )
    i_indices, j_indices = np.triu_indices(num_seqs, k=1)
    pos_distances = np.abs(mid_positions[i_indices] - mid_positions[j_indices])

    pos_distances[overlap_matrix[i_indices, j_indices]] = 0
    pos_distances[compensatory_matrix_dense[i_indices, j_indices]] = 0

    # fill distance matrix
    dist_matrix[i_indices, j_indices] = pos_distances
    # print("\ndist Matrix:")
    # print('\n'.join([' '.join(map(str, row)) for row in dist_matrix]))
    dist_matrix += dist_matrix.T  # Symmetrize the matrix
    print("finished calculating distance 4 step",datetime.now() )
    return dist_matrix, overlap_matrix, compensatory_matrix_dense


# recursive classification with top-down approach
def process_cluster_fully_overlap(root_clade, cluster_id, overlap_matrix, ids, fully_overlapping_clusters):
    stack = [(root_clade, cluster_id)] 
    while stack:
        current_clade, current_cluster_id = stack.pop()
        terminals = current_clade.get_terminals()
        member_ids = [terminal.name for terminal in terminals]
        # print("memeber_ids:",member_ids)
        member_ids = list(map(int, member_ids))
        # print("ids",ids)
        # print("memeber_ids:",member_ids)
        indices = [ids.index(name) for name in member_ids]
        num_sequences = len(indices)
        if num_sequences > 1:
            pair_indices = [(i, j) for idx_i, i in enumerate(indices)
                            for idx_j, j in enumerate(indices) if idx_i < idx_j]
            # 1st situation: every sequence overlaps with all other sequences
            all_overlaps = all(overlap_matrix[i, j] for i, j in pair_indices)

        if all_overlaps:
            fully_overlapping_clusters.append((current_cluster_id, member_ids, indices))
            continue

        # add child clades to the stack for further processing
        for child_clade in current_clade.clades:  # Use `clades` to access direct children
            if not child_clade.is_terminal():
                stack.append((child_clade, current_cluster_id + 1))

def process_cluster_all_compensatory(root_clade, cluster_id, compensatory_matrix, ids, fully_compensatory_clusters):
    stack = [(root_clade, cluster_id)] 
    while stack:
        current_clade, current_cluster_id = stack.pop()
        terminals = current_clade.get_terminals()
        member_ids = [terminal.name for terminal in terminals]
        member_ids = list(map(int, member_ids))
        indices = [ids.index(name) for name in member_ids]
        num_sequences = len(indices)
        if num_sequences > 1:
            all_compensatory = all(
                any(compensatory_matrix[i, j] for j in indices if i != j)
                for i in indices
            )

        if all_compensatory:
            fully_compensatory_clusters.append((current_cluster_id, member_ids, indices))
            continue

        # add child clades to the stack for further processing
        for child_clade in current_clade.clades:  # Use `clades` to access direct children
            if not child_clade.is_terminal():
                stack.append((child_clade, current_cluster_id + 1))

# recursive classification with top-down approach
def process_cluster_any_overlap(root_clade, cluster_id, overlap_matrix, ids, partially_overlapping_clusters):
    stack = [(root_clade, cluster_id)] 
    while stack:
        current_clade, current_cluster_id = stack.pop()
        terminals = current_clade.get_terminals()
        member_ids = [terminal.name for terminal in terminals]
        member_ids = list(map(int, member_ids))
        indices = [ids.index(name) for name in member_ids]
        num_sequences = len(indices)
        if num_sequences > 1:
            any_overlap = all(
                any(overlap_matrix[i, j] for j in indices if i != j)
                for i in indices
            )
        else:
            any_overlap = False
        # if parent_compensatory:
        #     # Skip processing below clusters 
        #     return

        if any_overlap:
            partially_overlapping_clusters.append((current_cluster_id, member_ids, indices))
            continue

        # add child clades to the stack for further processing
        for child_clade in current_clade.clades:  
            if not child_clade.is_terminal():
                stack.append((child_clade, current_cluster_id + 1))

def process_cluster_singleton(root_clade, cluster_id, ids, singularity_cluster):
    logging.info("start singleton tree")
    stack = [(root_clade, cluster_id)] 

    while stack:
        current_clade, current_cluster_id = stack.pop()
        terminals = current_clade.get_terminals()
        member_ids = [terminal.name for terminal in terminals]
        member_ids = list(map(int, member_ids))
        indices = [ids.index(name) for name in member_ids]
        num_sequences = len(indices)

        single = num_sequences == 1

        if single:
            singularity_cluster.append((current_cluster_id, member_ids, indices))
            continue 

        for child_clade in current_clade.clades: 
            if not child_clade.is_terminal():
                stack.append((child_clade, current_cluster_id + 1))

# merge overlapping x-axis intervals using IntervalTree
def merge_intervals(intervals):
    tree = IntervalTree(Interval(start, end) for start, end in intervals)
    tree.merge_overlaps()
    return sorted((iv.begin, iv.end) for iv in tree)

# function to plot each dataset
def plot_dataset(dataset, dataset_idx, output_folder):
    # convert each tuple (string values) into a list (integer values)
    chrom = dataset[0][0]
    dataset = [[int(start), int(end), int(line_count)] for chrom, start, end, line_count in dataset]
    plt.figure(figsize=(12, 4))
    x_min = min(int(item[0]) for item in dataset)
    x_max = max(int(item[1]) for item in dataset)
    parameter = 0.05
    # create unique x-limits based on dataset
    important_ranges = [(item[0], item[1]) for item in dataset]
    # logging.info("dataset_idx")
    # logging.info(dataset_idx)
    # logging.info("draw data start")
    # logging.info(dataset)
    # logging.info(important_ranges)
    merged_xlims = merge_intervals(important_ranges)
    # logging.info(merged_xlims)
    # logging.info("draw data end")
    distance = 0
    for start, end, line_count in dataset:
        if merged_xlims and int(start) == merged_xlims[0][0]:
            # logging.info("start")
            # logging.info(start)
            distance = 0
            merged_xlims.pop(0)
        y_positions = [distance + parameter * (i + 1) for i in range(int(line_count))]
        # print(y_positions)
        # random_color = np.random.rand(3,)
        for y in y_positions:
            plt.hlines(y, start, end, colors="black", linestyles='solid', linewidth=1)
        distance = y_positions[-1]
    plt.axis("off")
    plt.tight_layout()

    #save the plot to a buffer as a grayscale image
    buffer_path = f"{output_folder}/temp_plot_{dataset_idx}_{chrom}_{x_min}_{x_max}.png"
    plt.savefig(buffer_path, format='png', dpi=300, bbox_inches='tight')
    plt.close()

    #open the saved image and post-process it to binary (0 or 255)
    img = Image.open(buffer_path).convert('L')
    binary_img = img.point(lambda x: 255 if x > 0 else 0, mode='1')

    #Save the binary image
    binary_img.save(f"{output_folder}/binary_plot_{dataset_idx}_{chrom}_{x_min}_{x_max}.png")
    ##new code end

    # set labels and title for each figure
    # plt.xlabel("Genomic Positions")
    # plt.ylabel(f"Regions")
    # plt.title(f"Genomic Regions with Multiple Lines (Dataset {dataset_idx + 1})")
    # plt.grid(axis='x', linestyle='--', alpha=0.7)
    # plt.tight_layout()
    # Save the figure with the specified filename
    # filename = os.path.join(output_folder, f'middle_results/compensatory_cluster_id_{dataset_idx + 1}.png')
    # # filename = f'cluster_id_{dataset_idx + 1}.png'
    # plt.savefig(filename)
    # plt.close() 
    # plt.show()

def extract_features(img_path, model):
    img = image.load_img(img_path, target_size=(224, 224))
    img_data = image.img_to_array(img)
    print(img_data)
    img_data = np.expand_dims(img_data, axis=0)
    img_data = preprocess_input(img_data)
    features = model.predict(img_data)
    return features.flatten()

def unsupervise_learning(image_paths, output_folder):
    # use pretrained vgg16 models
    base_model = VGG16(weights='imagenet', include_top=False)
    print("image_paths",image_paths)
    binary_plot_files = [img_path for img_path in listdir(image_paths) if "binary_plot" in img_path]
    print("binary_plot_files:",binary_plot_files)
    features = [extract_features(path.join(image_paths, img_path), base_model) for img_path in binary_plot_files]
    print("features",features)
    features = np.array(features)
    print(f"Number of features extracted: {len(features)}")
    print(f"Features shape: {np.array(features).shape}")
    if len(features) > 0:
        print(f"First feature shape: {features[0].shape}")

    # pca lower the dimension
    n_components = min(len(features),len(binary_plot_files)) if min(len(features),len(binary_plot_files)) < 5000 else 5000
    pca = PCA(n_components=n_components)
    features_reduced = pca.fit_transform(features)

    # perform K-Means clustering
    ##first, calculate the optimal n_clusters
    silhouette_scores = []
    max_clusters = min(30,len(binary_plot_files)-1) 
    # calculate silhouette scores for different cluster sizes
    for n in range(2, max_clusters + 1):
        kmeans = KMeans(n_clusters=n, random_state=42, n_init='auto')
        kmeans.fit(features_reduced)
        score = silhouette_score(features_reduced, kmeans.labels_)
        silhouette_scores.append(score)
    optimal_clusters = 2 + np.argmax(silhouette_scores)
    print(f"Optimal number of clusters: {optimal_clusters}")

    kmeans = KMeans(n_clusters=optimal_clusters, random_state=42, n_init='auto')
    kmeans.fit(features_reduced)
    # get the cluster centroids
    centroids = kmeans.cluster_centers_
    print("centroid", centroids.shape)
    # get cluster centers and labels
    labels = kmeans.labels_

    # Step 3: visualize the results
    plt.figure(figsize=(8, 6))
    plt.scatter(features_reduced[:, 0], features_reduced[:, 1], c=labels, cmap='viridis', s=50, alpha=0.6, label='Data Points')
    plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=200, marker='X', label='Cluster Centers')
    plt.title('KMeans Clustering')
    plt.xlabel('Feature 1')
    plt.ylabel('Feature 2')
    plt.savefig(f"{output_folder}/KMeans_plot.png", format='png', dpi=300, bbox_inches='tight')
    plt.close()
    # calculate distances from each image to the centroids
    distances = euclidean_distances(features_reduced, centroids)
    # find the closest image for each centroid
    closest_images = np.argmin(distances, axis=0)
    print("closest_images",closest_images)

    return closest_images
   
def convert_labels_to_strings(tree):
    #ensure all node labels in the tree are strings
    for node in tree.preorder():
        if node.name is not None:
            node.name = str(node.name)
    return tree

def process_clusters(cluster_type, clusters, dataset_key, output_file_path, list_dict, bedfile_partial_df, sliding_window_id):
    """
    Processes different types of clusters and writes the results to a file.

    Parameters:
    - cluster_type: A string indicating the cluster type (for logging/debugging).
    - clusters: The cluster data containing (cluster_id, member_ids, indices).
    - dataset_key: The key in list_dict to store extracted data.
    - output_file_path: The file path to write output.
    - list_dict: Dictionary containing various datasets and file paths.
    - bedfile_partial_df: DataFrame containing BED file data.
    - sliding_window_id: Identifier for the sliding window.
    """
    with open(output_file_path, 'a') as output_file:
        for cluster_id, member_ids, indices in clusters:
            selected_bed_lines = get_bed_lines_by_indices(bedfile_partial_df, indices)
            selected_lines = selected_bed_lines.apply(lambda row: '\t'.join(map(str, row)), axis=1).tolist()
            
            # extract relevant columns from each line
            data = [(line.split('\t')[0], line.split('\t')[1], line.split('\t')[2], line.split('\t')[6]) for line in selected_lines]
            if data:
                list_dict[dataset_key].append(data)

            # format output text
            bed_lines_str = '\n'.join(selected_lines)
            output_file.write(f"\ncluster id:{sliding_window_id}_{cluster_id}, in this cluster's all bed ids : {member_ids}")
            output_file.write(f"\ncluster id:{sliding_window_id}_{cluster_id}, in this cluster's all bed lines :\n{bed_lines_str}\n")



def analyze_tree(bedfile_partial_df, ids, dist_matrix, method, overlap_matrix, compensatory_matrix, output_folder, list_dict, sliding_window_id):
    # def analyze_tree(output_folder):
    condensed_dist_matrix = squareform(dist_matrix)
    if condensed_dist_matrix.size == 0:
        return
    Z = linkage(condensed_dist_matrix, method=method)
    #save tree
    tree = TreeNode.from_linkage_matrix(Z, id_list=ids)
    tree = convert_labels_to_strings(tree)
    tree_path = os.path.join(output_folder, "tree.newick")
    # # check if the file exists and delete it
    # if os.path.exists(tree_path):
    #     os.remove(tree_path)
    #     print(f"Deleted existing file: {tree_path}")
    with open(tree_path, 'w') as f:
        tree.write(f, format='newick')
    print("finished building tree", datetime.now())

    # # load tree
    tree = Phylo.read(tree_path, 'newick')
    # classify clusters based on overlap and sequence symmetry
    fully_overlapping_clusters = []
    fully_compensatory_clusters = []
    partially_overlapping_clusters = []
    singularity_cluster = []
    # other_situation = []
    ids_int = list(map(int, ids))
    if overlap_matrix.size >=2 :
        process_cluster_fully_overlap(tree.root, 0, overlap_matrix, ids_int, fully_overlapping_clusters)
        process_cluster_any_overlap(tree.root, 0, overlap_matrix, ids_int, partially_overlapping_clusters)
    if compensatory_matrix.size >=2 :
        process_cluster_all_compensatory(tree.root, 0, compensatory_matrix, ids_int, fully_compensatory_clusters)
    process_cluster_singleton(tree.root, 0, ids_int, singularity_cluster)
    print("finished calculating tree 3 step",datetime.now() ) 
    logging.info("finished calculating tree 3 step")

    # Call the function for each cluster type
    process_clusters("Full Overlap", fully_overlapping_clusters, "dataset_full_overlap", list_dict["output_file_full_overlap"], list_dict, bedfile_partial_df, sliding_window_id)
    process_clusters("Compensatory", fully_compensatory_clusters, "dataset_compensatory", list_dict["output_file_compensatory"], list_dict, bedfile_partial_df, sliding_window_id)
    process_clusters("Partial Overlap", partially_overlapping_clusters, "dataset_partial_overlap", list_dict["output_file_partial_overlap"], list_dict, bedfile_partial_df, sliding_window_id)
    process_clusters("Singleton", singularity_cluster, "dataset_singleton", list_dict["output_file_singleton"], list_dict, bedfile_partial_df, sliding_window_id)
    # # Convert lists to sets and find the intersection
    # qseqid_list.sort()
    # intersection = set(mirdeep2_qseqid_list) & set(qseqid_list)
    # # Calculate the overlap rate
    # overlap_rate = 0
    # if len(mirdeep2_qseqid_list) > 0:
    #     overlap_rate = len(intersection) / len(mirdeep2_qseqid_list) # Relative to the smaller list
    # else:
    #     print("len(mirdeep2_qseqid_list) = 0")
    #     overlap_rate = -1
    # overlap_count = len(intersection)
    # output_file2.write(f"\nmirdeep2_qseqid_list: {mirdeep2_qseqid_list}\n qseqid_list:{qseqid_list}\n")
    # output_file2.write(f"\nmirdeep2_qseqid_list count: {len(mirdeep2_qseqid_list)}\n qseqid_list count:{len(qseqid_list)}\n")
    # output_file2.write(f"\nmirdeep2 qseqid including count: {overlap_count:.2f}, including qseqid: {intersection},mirdeep2 qseqid including rate: {overlap_rate * 100:.2f}%\n")
    # # # Plot each dataset in separate figures


def final_step(list_dict, output_folder):
    tmp_output_folder = "/tmp/mzhou10/image_results"
    if os.path.exists(tmp_output_folder):
        shutil.rmtree(tmp_output_folder)  # Remove all contents of the folder
    os.makedirs(tmp_output_folder) 
    # Check if the folder is readable
    if os.access(tmp_output_folder, os.R_OK):
        print(f"At the beginning, the folder {tmp_output_folder} is readable.")
    else:
        print(f"At the beginning, the folder {tmp_output_folder} is not readable.")
        sys.exit(1)
    
    combined_dataset = list_dict["dataset_full_overlap"] + list_dict["dataset_compensatory"] + list_dict["dataset_partial_overlap"] + list_dict["dataset_singleton"]
    if combined_dataset:
        for idx, dataset in enumerate(combined_dataset):
            plot_dataset(dataset, idx, tmp_output_folder)
        # image_path = os.path.join(tmp_output_folder, "middle_results")
        closest_images = unsupervise_learning(tmp_output_folder, output_folder)
        # Write centroid-to-image mapping to the output file
        output_file = os.path.join(output_folder, 'image_cluster_mapping.txt')
        binary_plot_files = [img_path for img_path in listdir(tmp_output_folder) if "binary_plot" in img_path]
        image_names =[]
        with open(output_file, 'w') as f:
            for cluster_idx, img_idx in enumerate(closest_images):
                img_name_1 = os.path.basename(binary_plot_files[img_idx])  # Get the image name
                img_name_2 = img_name_1.replace("binary", "temp") 
                image_names.append(img_name_1)
                image_names.append(img_name_2)
                f.write(f"Cluster {cluster_idx}: Closest image is {img_name_1},{img_name_2}\n")
        print(f"Cluster-to-image mapping written to {output_file}")

        for file_name in image_names:
            source_file = os.path.join(tmp_output_folder, file_name)
            if os.path.isfile(source_file):
                shutil.copy(source_file, path.join(output_folder,"middle_results"))
                print(f"Copied: {file_name}")
            else:
                print(f"File not found: {file_name}")

    # remove the original folder
    shutil.rmtree(tmp_output_folder)
    # Check if the folder is readable after removal
    if os.access(tmp_output_folder, os.R_OK):
        print(f"After finishing, the folder {tmp_output_folder} is readable.")
    else:
        print(f"After finishing, the folder {tmp_output_folder} is not readable.")
        sys.exit(1)


def classify_pattern(bed_file, output_folder):
    # print("1")
    log_filename = os.path.join(output_folder, "my_log_file.log")  # Specify the log file
    logging.basicConfig(
        filename=log_filename,  # Log to this file
        level=logging.INFO,  # Set log level
        format="%(asctime)s - %(levelname)s - %(message)s"  # Define log format
    )
    # # Example usage
    # positions = [
    #     [1, 5],  # Sequence 1
    #     [6, 10], # Sequence 2
    #     [11, 15] # Sequence 3
    # ]
    # gap_matrix = calculate_gap_matrix(positions)
    # print("Gap Matrix:")
    # print(gap_matrix)
    ##! need add slide window way to read bed file
    # write the classification results
    list_names = [ "dataset_full_overlap", "dataset_compensatory", "dataset_partial_overlap","dataset_singleton",  
                   "output_file_full_overlap", "output_file_compensatory", "output_file_partial_overlap", "output_file_singleton" ]
    list_dict = {name: [] if 'dataset' in name else "" for name in list_names}
    list_dict["output_file_full_overlap"] = os.path.join(output_folder, f"fully_overlapping_clusters.txt")
    list_dict["output_file_compensatory"] = os.path.join(output_folder, f"fully_compensatory_clusters.txt")
    list_dict["output_file_partial_overlap"] = os.path.join(output_folder, f"partially_overlapping_clusters.txt")
    list_dict["output_file_singleton"] = os.path.join(output_folder, f"singularity_cluster.txt")
    # list_dict["output_file_path5"] = os.path.join(output_folder, "other_situation.txt")
    sliding_window_id = 0

    for ids, sequences, positions, bedfile_partial_df in process_bed_file_in_ranges(bed_file):
        # ids, sequences, positions = read_bed_file(bed_file)
        # print("1")
        # logging.info("1")
        logging.info(bedfile_partial_df)
        columns_to_show = ["start", "end"]  # Replace with the columns you want
        logging.info(bedfile_partial_df.loc[:, columns_to_show])
        if ids and sequences and positions and not bedfile_partial_df.empty:
            dist_matrix, overlap_matrix, compensatory_matrix = compute_distance_matrix(sequences, positions, output_folder)
            print("2")
            # print("\ndist Matrix now:")
            # print('\n'.join([' '.join(map(str, row)) for row in dist_matrix]))
            # print("\ndist matrix size now",dist_matrix.size)
            if dist_matrix.size != 0:
                analyze_tree(bedfile_partial_df, ids, dist_matrix, 'average', overlap_matrix, compensatory_matrix, output_folder, list_dict, sliding_window_id)
                # print("3")
                print("dist_matrix not empty")
                sliding_window_id += 1
            else:
                print("dist_matrix empty")
                print(ids, sequences, positions, bedfile_partial_df, dist_matrix)  
    
    final_step(list_dict, output_folder)
    print("4")

def parse_arguments():
    parser = argparse.ArgumentParser(description="claassify mapping patterns")
    parser.add_argument('-input_bedfile', required=False,
                        type=str, help='input bed file in additional step')
    parser.add_argument('-output_folder', required=True,
                        type=str, help='output resutls folder path')

    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input_bedfile and args.output_folder:
        time_start_s = datetime.now()
        classify_pattern(args.input_bedfile, args.output_folder)
        time_end_s = datetime.now()
        time_c = time_end_s - time_start_s
        print('time cost', time_c, 's')


if __name__ == "__main__":
    main()