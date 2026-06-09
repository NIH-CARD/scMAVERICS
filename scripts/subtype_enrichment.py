import scanpy as sc
import decoupler as dc
import pandas as pd

# Import adata object
adata = sc.read_h5ad(work_dir + '/atlas/07_polished_anndata_rna.h5ad')

# Import subtype specific gene sets
subtype_df = pd.read_csv(work_dir + '/input/MG_subtype_genes.csv')

# Define the celltype and extract
celltype = 'MG'
subtype_adata = adata[adata.obs['celltype'] == celltype]

# Compute per-cell gene set enrichment
dc.mt.ulm(
    data = MG_adata,
    net = MG_subtype_df,
    raw = False,
    layer = 'log-norm', 
    tmin=2
)

# Extract the AnnData object with per-subtype scores
score = dc.pp.get_obsm(subtype_adata, key="score_ulm")

# Plot score per cell
sc.pl.umap(
    score, 
    color=subtype_df['source'].drop_duplicates().to_list + ["leiden_2"], 
    cmap="Reds",
    ncols = 3,
    show = False)
plt.savefig(work_dir + f'/figures/{work_dir}_subtype_prediction_scores.svg')

# Compute leiden clustering 
sc.tl.leiden(MG_adata, resolution=1, key_added = 'leiden', flavor = 'igraph')

# For each leiden cluster, calculate the most enriched subtype score
df = dc.tl.rankby_group(adata=score, groupby="leiden", reference="rest", method="t-test_overestim_var")
df = df[df["stat"] > 0]

# Assign subtype score
dict_ann = df[df["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
subtype_adata.obs["subtype"] = [dict_ann[x] for x in subtype_adata.obs["leiden"]]

# Write out AnnData object
subtype_adata.write_h5ad(work_dir + '/data/celltype/{celltype}_subtype.h5ad', compression = 'gzip')
