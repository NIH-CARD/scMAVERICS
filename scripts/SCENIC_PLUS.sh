#!/bin/bash

#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --partition=multinode
#SBATCH --time 48:00:00
#SBATCH --output=logs/SCENICPLUS-%j.out

module load bedtools/2.30.0 

REGION_BED="/data/CARD_singlecell/SN_atlas/data/pycisTopic/consensus_regions.bed"
GENOME_FASTA="/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
CHROMSIZES="/data/CARD_singlecell/SN_atlas/data/SCENICPLUS/hg38.chrom.sizes"
DATABASE_PREFIX="SN_atlas"
SCRIPT_DIR="/data/CARD_singlecell/SN_atlas/data/SCENICPLUS/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38.10x_brain.fa \
        1000 \
        yes

OUT_DIR="/data/CARD_singlecell/SN_atlas/data/SCENICPLUS"
CBDIR="${OUT_DIR}/aertslab_motif_collection/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/hg38.10x_brain.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 64
