import anndata as ad
import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import sys

# Import AnnData files
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)
filtered_adata = sc.read_h5ad(snakemake.input.hvg_rna_anndata)

# Save the .obsm['X_scvi']
adata.obsm['X_scvi'] = filtered_adata.obsm['X_scvi']

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
adata.obs['atlas_identifier'] = adata.obs.index.to_list()

# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata, min_dist=0.3)

# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(adata, resolution=2, key_added='leiden_2')
sc.tl.leiden(adata, key_added='leiden')
sc.tl.leiden(adata, resolution=.5, key_added='leiden_05')

# Save the anndata object
adata.write_h5ad(snakemake.output.merged_rna_anndata, compression='gzip')