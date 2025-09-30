#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00

python /data/CARD_singlecell/SN_atlas/scripts/rsid_linkage_disequilibrium.py