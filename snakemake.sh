#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:2
#SBATCH --array=0-6

module purge
module load apptainer
module load snakemake/7.7.0

# Pull profile, this will only run once, and is required for running on Biowulf
git clone https://github.com/NIH-HPC/snakemake_profile.git

# Pull the containers
apptainer pull envs/single_cell_cpu.sif oras://quay.io/adamcatchingdti/single_cell_cpu:0.4
apptainer pull envs/single_cell_gpu.sif oras://quay.io/adamcatchingdti/single_cell_gpu:0.8

apptainer pull envs/scenicplus.sif docker://litd/docker-scenicplus:latest
apptainer pull envs/decoupler.sif docker://deeenes/omnipath-decoupler

# Load singularity
module load singularity/4.1.5

# Bind external directories on Biowulf
. /usr/local/current/singularity/app_conf/sing_binds

# RUN SCRIPT
snakemake --cores all --profile snakemake_profile --use-singularity 
