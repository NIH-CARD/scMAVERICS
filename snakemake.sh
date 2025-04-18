#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 96:00:00

# SETUP
module purge
module load apptainer
module load snakemake/7.7.0
module load singularity

# Pull profile, this will only run once, and is required for running on Biowulf
git clone https://github.com/NIH-HPC/snakemake_profile.git

# Pull the containers
# - Public 
apptainer pull envs/snapatac2.sif oras://quay.io/adamcatchingdti/snapatac2
apptainer pull envs/single_cell_gpu.sif oras://quay.io/adamcatchingdti/single_cell_gpu:0.8
apptainer pull envs/decoupler.sif oras://quay.io/adamcatchingdti/decoupler.sif:0.8
# - Personal
apptainer pull envs/scenicplus.sif docker://litd/docker-scenicplus:latest

# Bind external directories on Biowulf
. /usr/local/current/singularity/app_conf/sing_binds

# Update permissions on the bash scripts 
chmod 777 scripts/rna_model.sh
chmod 777 scripts/cellbender_array.sh
chmod 777 scripts/atac_model.sh

# RUN SCRIPT
snakemake \
    --cores all \
    --profile snakemake_profile \
    --use-singularity
    # --use-singularity annotate 
    #-n
    #--unlock