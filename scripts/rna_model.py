import anndata as ad
import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import sys

print(torch.cuda.is_available())

scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')

print(sys.argv)

# Read in AnnData atlas object
filtered_adata = ad.read_h5ad(sys.argv[1])

# Setup SCVI on the data layer
scvi.model.SCVI.setup_anndata(
    filtered_adata, layer="log-norm", batch_key=sys.argv[2])

# Add the parameters of the model
model = scvi.model.SCVI(
    filtered_adata, 
    dispersion="gene-batch", 
    n_layers=2, 
    n_latent=30, 
    gene_likelihood="nb"
)

# Train the model
model.train(
    max_epochs=1000,  
    early_stopping=True,
    early_stopping_patience=20
)

model.save(sys.argv[5], overwrite=True)

# Save the anndata object
filtered_adata.write_h5ad(sys.argv[4], compression='gzip')

"""# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(sys.argv[3], index=False)

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
adata.obs['atlas_identifier'] = adata.obs.index.to_list()
adata.obsm['X_scvi'] = model.get_latent_representation()

# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata, min_dist=0.3)
# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(adata, resolution=2, key_added='leiden_2')
sc.tl.leiden(adata, key_added='leiden')
sc.tl.leiden(adata, resolution=.5, key_added='leiden_05')

# Save the anndata object
adata.write_h5ad(sys.argv[4], compression='gzip')

model.save(sys.argv[5], overwrite=True)"""
