import scanpy as sc
import pandas as pd
import numpy as np
import sys
import muon as mu

# Import Muon files
mdata = mu.read(snakemake.input.multiome_object)
filtered_mdata = mu.read(snakemake.input.hvg_multiome_anndata)

# Save the .obsm['X_scvi']
mdata.obsm['X_multivi'] = filtered_mdata.obsm['X_multivi']

# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
sc.pp.neighbors(mdata, use_rep='X_multivi')
sc.tl.umap(mdata, min_dist=0.3)
# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(mdata, key_added='leiden', flavor = 'igraph')

# Save the anndata object
mdata.write(snakemake.output.multiome_object, compression='gzip')