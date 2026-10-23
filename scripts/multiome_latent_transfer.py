import scanpy as sc
import pandas as pd
import numpy as np
import sys
import muon as mu
import snapatac2 as snap

# Import Muon files
mdata = mu.read(snakemake.input.multiome_object)
filtered_mdata = mu.read(snakemake.input.hvg_multiome_anndata)

# Save the .obsm['X_scvi']
mdata.obsm['X_multivi'] = filtered_mdata.obsm['X_multivi']
mdata.obsm['X_umap'] = filtered_mdata.obsm['X_umap']
mdata.obsm['X_multivi'] = filtered_mdata.obsm['X_multivi']

mdata.obsp['connectivities'] = filtered_mdata.obsp['connectivities']
mdata.obsp['distances'] = filtered_mdata.obsp['distances']


mdata.obs['leiden_2'] = filtered_mdata.obs['leiden_2']
mdata.obs['leiden'] = filtered_mdata.obs['leiden']


# Port leiden clusters to rna and atac modalities
barcode2leiden = mdata.obs['leiden'].to_dict()
mdata.mod['rna'].obs['leiden'] = [barcode2leiden[x] for x in mdata.mod['rna'].obs_names]
mdata.mod['atac'].obs['leiden'] = [barcode2leiden[x] for x in mdata.mod['atac'].obs_names]

"""# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
print('running nearest neighbors')
sc.pp.neighbors(mdata, use_rep='X_multivi')
print('running umap')
snap.tl.umap(mdata, use_rep='X_multivi')
# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
print('running leiden')
sc.tl.leiden(mdata, key_added='leiden', flavor = 'igraph')
sc.tl.leiden(mdata, key_added='leiden_2', resolution = 2, flavor = 'igraph')"""

mdata = mdata[~(mdata.obs['leiden'].isin(['10', '17', '29', '31', '34', '35']))].copy()
mdata

# Save the anndata object
mdata.write(snakemake.output.multiome_object, compression='gzip')