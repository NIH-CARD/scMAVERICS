#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1

"""
This takes all of the transferred Cellranger output files and runs cellbender
"""

# Load modules
module load cellbender/0.3.2
module load CUDA/12.1

# The temp directory is hardcoded in CellBender for where checkpoint files are output
export TMPDIR=$2
echo $TMPDIR

# Iterate through the array of sample directories
cd $TMPDIR; pwd; cellbender remove-background --input $1 --output $3 --checkpoint "$TMPDIR"/ckpt.tar.gz --cuda; cd ..