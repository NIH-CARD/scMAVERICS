#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=8
#SBATCH --time 24:00:00
#SBATCH --gres=lscratch:2000
#SBATCH --output=/data/CARD_singlecell/SN_atlas/logs/homer.out


module load homer

annotatePeaks.pl data/pycisTopic/consensus_regions.bed hg38 > data/homer/annotated_consensus_regions.txt