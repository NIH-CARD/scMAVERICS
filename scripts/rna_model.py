import anndata as ad
import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import sys

print(torch.cuda.is_available())

scvi.settings.seed = sys.argv[6]
torch.set_float32_matmul_precision('high')

# Read in AnnData atlas object
adata = ad.read_h5ad(sys.argv[1])

# Setup SCVI on the data layer
scvi.model.SCVI.setup_anndata(
    adata, layer="counts", batch_key=sys.argv[2])

# Add the parameters of the model
model = scvi.model.SCVI(
    adata, 
    dispersion="gene-batch", 
    n_layers=sys.argv[7], 
    n_latent=sys.argv[8], 
    gene_likelihood="nb"
)

# Train the model
model.train(
    max_epochs=sys.argv[9],
    accelerator=sys.argv[10],  
    early_stopping=True,
    early_stopping_patience=20
)

# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(sys.argv[3], index=False)

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
adata.obsm['X_scvi'] = model.get_latent_representation()

# Save the anndata object
adata.write_h5ad(sys.argv[4], compression='gzip')

model.save(sys.argv[5], overwrite=True)
