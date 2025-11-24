# Load in modules
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import pickle



# Import ATAC data
atac = sc.read_h5ad(snakemake.input.atac_anndata)

# Filter for only control data
control_atac = atac[atac.obs[snakemake.params.disease_param] == snakemake.params.control]

# Get marker peaks
marker_peaks = snap.tl.marker_regions(
    data = control_atac,
    groupby = snakemake.params.cell_type, 
    pvalue = 0.01
)

# Read in TF motifs (must be MEME format)
TF_motifs = snap._snapatac2.read_motifs(snakemake.input.TF_motifs)

# Motif enrichment, use CIS-BP TF list
motifs = snap.tl.motif_enrichment(
    motifs = jaspar_2024,
    regions = marker_peaks,
    genome_fasta = snakemake.input.ref_genome
)

# Read each cell type into a motif
motif_df = pd.DataFrame()
for key in motifs.keys():
    cell_motif = motifs[key].to_pandas()
    cell_motif[snakemake.params.cell_type] = key
    motif_df = pd.concat([motif_df, cell_motif])

# Export motif enrichment DataFrame
motif_df.to_csv(snakemake.output.motif_enrichment)