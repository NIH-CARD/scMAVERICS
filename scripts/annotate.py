import numpy as np
import scanpy as sc
import decoupler as dc
import scvi
import pandas as pd

# Open the RNA merged and filtered
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Create the DataFrame of canonical gene markers (This can be expanded)
marker_genes = pd.read_csv(snakemake.input.gene_markers)
marker_gene_df = pd.DataFrame(marker_genes)

# Run over-represenation analysis based on cell markers
# provided in the marker_gene_df DataFrame.
dc.run_ora(
    mat=adata,
    net=marker_gene_df,
    source='cell type',
    target='official gene symbol',
    min_n=1,
    verbose=True,
    use_raw=False
)

# Create a mini AnnData object with the over-represenation
# analysis estimate (p-value of given cell marker)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

# Convert the ORA AnnData object to numpy array to rank
# which cell type for each leiden cluster
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e
df = dc.rank_sources_groups(
    acts, 
    groupby='leiden_2', 
    reference='rest', 
    method='t-test_overestim_var'
    )

# Apply the best ranked cell type to a cluster-celltype dictionary
annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

# Apply the dictionary to the AnnData object
adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden_2']]

# Save the cell barcode, cluster, cell-type, and batch values to a .csv
adata.obs[['atlas_identifier', 'leiden_2', 'cell_type', 'Use_batch']].to_csv(snakemake.output.cell_annotate, index=False)

# Save the annotated AnnData object
adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip')