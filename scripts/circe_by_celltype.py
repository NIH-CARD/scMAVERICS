import numpy as np
import circe as ci
import scanpy as sc
import scipy as sp

adata = sc.read_h5ad(snakemake.input.celltype_atac)
adata.var['start'] = adata.var['start'].astype(int)
adata.var['end'] = adata.var['end'].astype(int)
adata.var['chr'] = adata.var['chromosome']

# Compute metacells
metacells = ci.metacells.compute_metacells(adata)

# Get co-accessibility scores
final_score = ci.sliding_graphical_lasso(
    adata,
    n_samples=50,
    n_samples_maxtry=100,
    max_alpha_iteration=500,
    verbose=True
)

adata.write_h5ad(snakemake.output.celltype_atac, compression='gzip')