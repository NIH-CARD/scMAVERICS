import anndata as ad
import muon as mu
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
mdata = mu.read(sys.argv[1])

# Setup SCVI on the data layer
scvi.model.MULTIVI.setup_mudata(
    mdata,
    modalities={
        "rna_layer": "rna",
        "atac_layer": "atac"
    }
)

# Train
mvi_model = scvi.model.MULTIVI(
    mdata,
    n_genes=mdata.mod["rna"].n_vars,
    n_regions=mdata.mod["atac"].n_vars,
)

mvi_model.train(
    accelerator="gpu",
    max_epochs=1000,
    lr=1e-3,
    early_stopping=True,
    batch_size=256,
    early_stopping_patience=20
)

# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(sys.argv[3], index=False)

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
adata.obsm['X_multivi'] = model.get_latent_representation()

# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
sc.pp.neighbors(adata, use_rep='X_multivi')
sc.tl.umap(adata, min_dist=0.3)
# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(adata, resolution=2, key_added='leiden_2')
sc.tl.leiden(adata, key_added='leiden')
sc.tl.leiden(adata, resolution=.5, key_added='leiden_05')

# Save the anndata object
adata.write_h5ad(sys.argv[4], compression='gzip')

model.save(sys.argv[5], overwrite=True)
