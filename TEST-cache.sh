#!/bin/bash
#SBATCH --gres=lscratch:100

# TEST
export TMPDIR=/lscratch/$SLURM_JOB_ID
echo "$SLURM_JOB_ID = {$SLURM_JOB_ID}"
echo "$TMPDIR = {$TMPDIR}"

# PROPOSAL
# replace this:
# export APPTAINER_CACHEDIR=/tmp/user/temporary-cache
# with this:
export APPTAINER_CACHEDIR=/lscratch/$SLURM_JOB_ID

# REFS:
# Apptainer: https://apptainer.org/docs/user/1.0/build_env.html
# Biowulf lscratch: https://hpc.nih.gov/docs/userguide.html#local