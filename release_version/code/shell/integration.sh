#!/bin/bash

module load blast/2.14
module load biodata/1.0

code_path=$1
mapping_species=$2
mapping_mirna=$3
output_folder=$4

###need modify python code
python3 $code_path/mirna_and_species_blast_overlap.py -species_folder $mapping_species -hairpin_folder $mapping_mirna -overlap_folder $output_folder

