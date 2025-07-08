#!/bin/bash
#SBATCH -o /data/CARD_singlecell/PFC_atlas/logs/01_atac_model.out
#SBATCH -e /data/CARD_singlecell/PFC_atlas/logs/01_atac_model.err
#SBATCH --partition=gpu 
#SBATCH --cpus-per-task=2
#SBATCH --mem=350000
#SBATCH --gres=gpu:v100x:2
#SBATCH --time=48:00:00

# This is so kludgy, but this the input files just need to be based through as arguments
input_file=/data/CARD_singlecell/PFC_atlas/atlas/03_merged_cistopic_atac.h5ad
sample_key='sample_id'
atac_model_history=/data/CARD_singlecell/PFC_atlas/data/model_elbo/atac_hmmr_model_history.csv
merged_atac_anndata=/data/CARD_singlecell/PFC_atlas/atlas/04_modeled_cistopic_atac.h5ad
atac_model=/data/CARD_singlecell/PFC_atlas/data/models/atac/

# Load module
module load singularity
# Run 
singularity run --nv --bind /data/CARD_singlecell/PFC_atlas envs/single_cell_gpu.sif python /data/CARD_singlecell/PFC_atlas/scripts/atac_model.py "${input_file}" "${sample_key}" "${atac_model_history}" "${merged_atac_anndata}" "${atac_model}"
