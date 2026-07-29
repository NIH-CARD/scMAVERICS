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
mdata.obsm['X_umap'] = filtered_mdata.obsm['X_umap']


mdata.obsp['distances']	= filtered_mdata.obsp['distances']
mdata.obsp['connectivities'] = filtered_mdata.obsp['connectivities']

mdata.uns['neighbors']	= filtered_mdata.uns['neighbors']
mdata.uns['umap'] = filtered_mdata.uns['umap']
mdata.uns['leiden']	= filtered_mdata.uns['leiden']

mdata.obs['leiden']	= filtered_mdata.obs['leiden']

# Save the anndata object
mdata.write(snakemake.output.multiome_object, compression='gzip')