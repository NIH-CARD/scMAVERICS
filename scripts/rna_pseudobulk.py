import scanpy as sc
import decoupler as dc

# Open rna
adata = sc.read_h5ad(snakemake.input.rna_anndata)

# Get pseudo-bulk profile, summing all raw counts by cell type and sample
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=snakemake.params.sample_key,
    groups_col=snakemake.params.separating_cluster,
    layer='counts',
    mode='sum'
)

# Filter out samples that have 
dc.pp.filter_samples(pdata, min_cells=snakemake.params.min_cells)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

dc.pp.filter_by_expr(
    adata=pdata,
    group="celltype",
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
)
dc.pp.filter_by_prop(
    adata=pdata,
    min_prop=0.1,
    min_smpls=2,
)

pdata.write_h5ad(snakemake.output.pseudo_rna, compression='gzip')