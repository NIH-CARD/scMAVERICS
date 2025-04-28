import anndata as ad
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import sys

# Read in AnnData atlas object
adata = ad.read_h5ad(snakemake.input.merged_rna_anndata)

# Select for the most variable genes
sc.pp.highly_variable_genes(
    adata, 
    layer='log-norm',
    n_top_genes=2000)

# Double check that no transcripts not found in cells are in the atlas
sc.pp.filter_genes(adata, min_cells=10)

# Define mitochondria and ribosome genes to remove
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Make a copy of the AnnData atlas that only contains variable genes
filtered_adata = adata[:, (adata.var['highly_variable']) & ~(adata.var['mt']) & ~(adata.var['rb'])].copy()

# Save the anndata object
filtered_adata.write_h5ad(snakemake.output.hvg_rna_anndata, compression='gzip')