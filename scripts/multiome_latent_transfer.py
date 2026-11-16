import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import sys
import muon as mu
import snapatac2 as snap
import os
# Force Numba and OpenMP to completely linear execution to stop the SIGSEGV
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMBA_NUM_THREADS"] = "1"


# Import Muon files
mdata = mu.read(snakemake.input.multiome_object)
filtered_mdata = mu.read(snakemake.input.hvg_multiome_anndata)

# Save the .obsm['X_scvi']
mdata.obsm['X_multivi'] = filtered_mdata.obsm['X_multivi']

joint_latent = mdata.obsm["X_multivi"] 

adata_sketch = ad.AnnData(obs=mdata.obs.copy())
adata_sketch.obsm["X_multivi"] = joint_latent

sc.tl.umap(
    adata_sketch, 
    init_pos="random",  
    min_dist=0.3,     # Tightens local packing to make clusters distinct
    maxiter=600,       # Gives the random initialization more cycles to resolve forces
    copy=False
)

mdata.obsm["X_umap"] = adata_sketch.obsm["X_umap"].copy()
mdata.uns["neighbors"] = adata_sketch.uns["neighbors"].copy()
# Sync modality configurations
mdata.update()

# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(mdata, key_added='leiden', flavor = 'igraph')
sc.tl.leiden(mdata, key_added='leiden_2', resolution = 2, flavor = 'igraph')

# Port leiden clusters to rna and atac modalities
barcode2leiden = mdata.obs['leiden'].to_dict()
mdata.mod['rna'].obs['leiden'] = [barcode2leiden[x] for x in mdata.mod['rna'].obs_names]
mdata.mod['atac'].obs['leiden'] = [barcode2leiden[x] for x in mdata.mod['atac'].obs_names]

# Save the anndata object
mdata.write(snakemake.output.multiome_object, compression='gzip')