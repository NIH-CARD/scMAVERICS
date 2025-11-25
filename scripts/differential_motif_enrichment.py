# Import modules
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap


# Import cell type atac anndata object
cell_type_atac = sc.read_h5ad(snakemake.input.cell_type_atac)

# Import cell type and disease specific DARs
cell_type_DARs = pd.read_csv(snakemake.input.output_DAR_data)

# Load in TF motifs
tf_motifs = snap._snapatac2.read_motifs(snakemake.input.TF_motifs)

# Return list of top and bottom 5% enriched peaks
DAR_counts = cell_type_DARs.shape[0] # Total number of peaks
top_05_num = int(DAR_counts * 0.05) # 5% of the total number

enriched_peaks = cell_type_DARs.sort_values('log2FoldChange', ascending=False).iloc[:top_05_num]['peak'].values
depleted_peaks = cell_type_DARs.sort_values('log2FoldChange', ascending=False).iloc[top_05_num:]['peak'].values

cell_disease_motifs = snap.tl.motif_enrichment(
    motifs=tf_motifs,
    regions={'up': enriched_peaks, 'down': depleted_peaks},
    background=cell_type_DARs['peak'].values,
    genome_fasta=ref_genome,
    method='hypergeometric'
)
# Split out up and down regulated motifs
up_df = cell_cre_disease_motifs['up'].to_pandas()
up_df['regulation'] = 'up'

down_df = cell_cre_disease_motifs['down'].to_pandas()
down_df = 'down'

cell_disease_motif_df = pd.concat([up_df, down_df])
cell_disease_motif_df['-log10(p-value)'] = -np.log10(cell_disease_motif_df['adjusted p-value'])

# Export dataframe
cell_disease_motif_df.to_csv(snakemake.output.differential_motif_dataframe, index=False)