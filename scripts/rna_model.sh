#!/bin/bash

# This is so kludgy, but this the input files just need to be based through as arguments
input_file=$1
sample_key=$2
model_history=$3
merged_rna_anndata=$4
model=$5

# Load module
module load singularity/4.1.5
# Run 
singularity run --nv --bind /data/CARD_singlecell/SN_atlas envs/single_cell_gpu.sif python /data/CARD_singlecell/SN_atlas/scripts/rna_model.py "${input_file}" "${sample_key}" "${model_history}" "${merged_rna_anndata}" "${model}"
