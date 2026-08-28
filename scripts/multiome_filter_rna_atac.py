import anndata as ad
import scanpy as sc
import pandas as pd

# Open the RNA merged and filtered
base_dir = '/data/CARD_singlecell/SN_atlas/'
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Add the consolidated cell-barcode 'atlas_identifier'
rna_samples = rna.obs['sample'].to_list()
rna_barcodes = rna.obs['cell_barcode'].to_list()
# Initialize cell-barcode 
rna_cell_barcode = []
for i in range(rna.n_obs):
    rna_cell_barcode.append(rna_samples[i] + '-' + rna_barcodes[i])

# Save the identifier
rna.obs['atlas_identifier'] = rna_cell_barcode

# Create the atlas identifier from the adataset obs
atac_df = atac.obs
atac_df['atlas_identifier'] = atac_df['index'] + '_' + atac_df['sample']
atac.obs = atac_df

# Subset RNA and ATAC objects based on the overlap of values
rna = rna[rna.obs['atlas_identifier'].isin(atac.obs['atlas_identifier'])].copy()
atac = atac[atac.obs['atlas_identifier'].isin(rna.obs['atlas_identifier'])].copy()

# Add the RNA observation data to the ATAC data
atac_df = pd.merge(
    left=atac_df,
    right=rna_df,
    left_on='atlas_identifier',
    right_on='atlas_identifier')
atac.obs = atac_df

# Write out the objects
rna.write_h5ad(snakemake.output.merged_rna_anndata, compression='gzip')
atac.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')