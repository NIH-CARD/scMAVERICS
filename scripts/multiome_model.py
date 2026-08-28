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

scvi.settings.seed = sys.argv[6]
torch.set_float32_matmul_precision('high')

# Read in AnnData atlas object
mdata = mu.read(sys.argv[1])
mdata
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
    accelerator=sys.argv[8],
    max_epochs=sys.argv[7],
    lr=1e-3,
    early_stopping=True,
    batch_size=256
)

# Extract the elbo plot of the model and save the values
elbo = mvi_model.history['elbo_train']
elbo['elbo_validation'] = mvi_model.history['elbo_validation']
elbo.to_csv(sys.argv[3], index=False)

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
mdata.obsm['X_multivi'] = mvi_model.get_latent_representation()

# Save the anndata object
mdata.write_h5ad(sys.argv[4], compression='gzip')

mvi_model.save(sys.argv[5], overwrite=True)
