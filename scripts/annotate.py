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

# Create the DataFrame of canonical gene markers (This can be expanded)
marker_gene_df = pd.read_csv(snakemake.input.gene_markers)

# Run over-represenation analysis based on cell markers
# provided in the marker_gene_df DataFrame.
dc.mt.ulm(
    data=adata, 
    net=marker_genes_df.rename(columns={'cell type' : 'source', 'official gene symbol': 'target'}), tmin=1)

# Create a mini AnnData object with the over-represenation
# analysis estimate (p-value of given cell marker)
score = dc.pp.get_obsm(adata, key="score_ulm")
df = dc.tl.rankby_group(adata=score, groupby="leiden", reference="rest", method="t-test_overestim_var")
df = df[df["stat"] > 0]

ctypes_dict = df.groupby("group").head(1).groupby("group")["name"].apply(lambda x: list(x)).to_dict()

dict_ann = df[df["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()

adata.obs["celltype"] = adata.obs["leiden"].cat.rename_categories(dict_ann)
# Apply the best ranked cell type to a cluster-celltype dictionary


# Save the cell barcode, cluster, cell-type, and batch values to a .csv
adata.obs[['atlas_identifier', 'leiden', 'celltype', snakemake.params.seq_batch_key]].to_csv(snakemake.output.cell_annotate, index=False)


# Save the annotated AnnData object
adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip')