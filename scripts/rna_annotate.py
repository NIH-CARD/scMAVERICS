import numpy as np
import scanpy as sc
import decoupler as dc
import scvi
import pandas as pd

# Open the RNA merged and filtered
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

doublet_clusters = []
for cluster in adata.obs[snakemake.params.leiden_cluster].drop_duplicates():
    #print(cluster, adata[adata.obs['leiden'] == cluster].obs['doublet_score'].mean(), adata[adata.obs['leiden'] == cluster].obs['doublet_score'].median())
    if adata[adata.obs[snakemake.params.leiden_cluster] == cluster].obs['doublet_score'].median() > .05:
        doublet_clusters.append(cluster)

adata = adata[~adata.obs[snakemake.params.leiden_cluster].isin(doublet_clusters)].copy()

mito_clusters = []
for cluster in adata.obs[snakemake.params.leiden_cluster].drop_duplicates():
    #print(cluster, adata[adata.obs['leiden'] == cluster].obs['doublet_score'].mean(), adata[adata.obs['leiden'] == cluster].obs['doublet_score'].median())
    if adata[adata.obs[snakemake.params.leiden_cluster] == cluster].obs['pct_counts_mt'].median() > 2:
        mito_clusters.append(cluster)
mito_clusters

# Create the DataFrame of canonical gene markers (This can be expanded)
marker_gene_df = pd.read_csv(snakemake.input.gene_markers)

# Run over-represenation analysis based on cell markers
# provided in the marker_gene_df DataFrame.
dc.mt.ulm(
    adata,
    net = marker_gene_df,
    tmin=1,
)

# Convert the ORA AnnData object to numpy array to rank
score = dc.pp.get_obsm(adata, key="score_ulm")
df = dc.tl.rankby_group(adata=score, groupby=snakemake.params.leiden_cluster, reference="rest", method="t-test_overestim_var")

# Apply the best ranked cell type to a cluster-celltype dictionary
annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

# Apply the dictionary to the AnnData object
adata.obs['celltype'] = [annotation_dict[clust] for clust in adata.obs[snakemake.params.leiden_cluster]]

# Save the annotated AnnData object
adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip')