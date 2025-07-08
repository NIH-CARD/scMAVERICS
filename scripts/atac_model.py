import anndata as ad
import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import sys

# To make sure GPU is accessible
print(torch.cuda.is_available())

# Set up parameters
scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')

# Read in AnnData atlas object
adata = ad.read_h5ad(sys.argv[1])

# Make matrix sparse
adata.X = scipy.sparse.csr_matrix(adata.X.astype(np.float64)[:])

print("# regions before filtering:", adata.shape[-1])

# Create object to filter without removing 
filtered_adata = adata.copy()
# Compute the threshold: 5% of the cells
min_cells = int(filtered_adata.shape[0] * 0.05)
# Filter genes 
sc.pp.filter_genes(filtered_adata, min_cells=min_cells)

print("# regions after filtering:", filtered_adata.shape[-1])

# Setup POISSONVI on the data layer
scvi.external.POISSONVI.setup_anndata(filtered_adata) 

# Add the parameters of the model
model = scvi.external.POISSONVI(filtered_adata)

# Train the model
model.train(
    max_epochs=1000,
    accelerator='gpu', 
    early_stopping=True,
    early_stopping_patience=20,
    )

# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(sys.argv[3], index=False)

# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
adata.obs['atlas_identifier'] = adata.obs.index.to_list()
adata.obsm['X_poissonvi'] = model.get_latent_representation()

# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
sc.pp.neighbors(adata, use_rep='X_poissonvi')
sc.tl.umap(adata, min_dist=0.3)

# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(adata, resolution=2, key_added='leiden_2')
sc.tl.leiden(adata, key_added='leiden')
sc.tl.leiden(adata, resolution=.5, key_added='leiden_05')

# Save the anndata object
adata.write_h5ad(sys.argv[4], compression='gzip')

# Save the atac model
model.save(sys.argv[5], overwrite=True)