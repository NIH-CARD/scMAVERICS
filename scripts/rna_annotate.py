import numpy as np
import scanpy as sc
import decoupler as dc
import scvi
import pandas as pd

# Open the RNA merged and filtered
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

doublet_clusters = []
for cluster in adata.obs['leiden'].drop_duplicates():
    #print(cluster, adata[adata.obs['leiden'] == cluster].obs['doublet_score'].mean(), adata[adata.obs['leiden'] == cluster].obs['doublet_score'].median())
    if adata[adata.obs['leiden'] == cluster].obs['doublet_score'].median() > .05:
        doublet_clusters.append(cluster)

adata = adata[~adata.obs['leiden'].isin(doublet_clusters)].copy()

mito_clusters = []
for cluster in adata.obs['leiden'].drop_duplicates():
    #print(cluster, adata[adata.obs['leiden'] == cluster].obs['doublet_score'].mean(), adata[adata.obs['leiden'] == cluster].obs['doublet_score'].median())
    if adata[adata.obs['leiden'] == cluster].obs['pct_counts_mt'].median() > 2:
        mito_clusters.append(cluster)
mito_clusters

# Create the DataFrame of canonical gene markers (This can be expanded)
marker_gene_df = pd.read_csv(snakemake.input.gene_markers)

# Run over-represenation analysis based on cell markers
# provided in the marker_gene_df DataFrame.
dc.mt.ulm(
    adata,
    net = gene_marker_df,
    min_n = 1,
    use_raw = False,
    layer = 'log-norm'
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
adata.obs['celltype'] = [annotation_dict[clust] for clust in adata.obs['leiden_2']]

# Save the annotated AnnData object
adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip')