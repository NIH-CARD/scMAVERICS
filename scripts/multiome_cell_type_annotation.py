import numpy as np
import pandas as pd

import scanpy as sc
import decoupler as dc
import muon as mu

mdata = mu.read_h5mu(snakemake.input.multiome_object)

# Transfer leiden
barcode2leiden = mdata.obs[snakemake.params.leiden_cluster].to_dict()
mdata.mod['rna'].obs[snakemake.params.leiden_cluster] = [barcode2leiden[x] for x in mdata.mod['rna'].obs_names]

doublet_clusters = []
for cluster in mdata.obs[snakemake.params.leiden_cluster].drop_duplicates():
    #print(cluster, adata[adata.obs['leiden'] == cluster].obs['doublet_score'].mean(), adata[adata.obs['leiden'] == cluster].obs['doublet_score'].median())
    if mdata[mdata.obs[snakemake.params.leiden_cluster] == cluster].obs['rna:doublet_score'].median() > .05:
        doublet_clusters.append(cluster)

mdata = mdata[~mdata.obs[snakemake.params.leiden_cluster].isin(doublet_clusters)].copy()

# Create the DataFrame of canonical gene markers (This can be expanded)
marker_gene_df = pd.read_csv(snakemake.input.gene_markers)

# Use only the rna modality
adata = mdata.mod['rna']

# Run over-represenation analysis based on cell markers
# provided in the marker_gene_df DataFrame
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

# Transfer to mdata
leiden2celltype = dict(zip(adata.obs[snakemake.params.leiden_cluster], adata.obs['celltype']))
mdata.obs['celltype'] = [leiden2celltype[x] for x in mdata.obs[snakemake.params.leiden_cluster]]

# Save the annotated AnnData object
mdata.write_h5mu(filename=snakemake.output.merged_rna_anndata, compression='gzip')