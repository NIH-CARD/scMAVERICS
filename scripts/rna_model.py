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

# Read in AnnData atlas object
filtered_adata = ad.read_h5ad(sys.argv[1])

# Setup SCVI on the data layer
scvi.model.SCVI.setup_anndata(
    filtered_adata, layer="counts", batch_key=sys.argv[2])

# Add the parameters of the model
model = scvi.model.SCVI(
    filtered_adata, 
    dispersion="gene-batch", 
    n_layers=2, 
    n_latent=30, 
    gene_likelihood="nb"
)

print('Starting modeling')
# Train the model
model.train(
    max_epochs=1000,
    early_stopping=True,
    accelerator='gpu',
    early_stopping_patience=40
)

print('Done modeling')

# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(sys.argv[3], index=False)

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
filtered_adata.obsm['X_scvi'] = model.get_latent_representation()

# Save the anndata object
filtered_adata.write_h5ad(sys.argv[4], compression='gzip')

# Done modeling
model.save(sys.argv[5], overwrite=True)
