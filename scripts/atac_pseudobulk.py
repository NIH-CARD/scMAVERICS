import scanpy as sc
import decoupler as dc

# Open atac
adata = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Get pseudo-bulk profile, summing all raw counts by cell type and sample
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=snakemake.params.sample_key,
    groups_col=snakemake.params.separating_cluster,
    mode='sum'
)

# Filter out samples that have 
dc.pp.filter_samples(pdata, min_cells=snakemake.params.min_cells)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

# Save pseudobulked atac data
pdata.write_h5ad(snakemake.output.pseudo_atac, compression='gzip')
