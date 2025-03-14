import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np

print(torch.cuda.is_available())

scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')

# Read in AnnData atlas object
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Double check that no transcripts not found in cells are in the atlas
sc.pp.filter_genes(adata, min_cells=1)

# Select for the most variable genes
sc.pp.highly_variable_genes(
    adata, 
    layer='data',
    n_top_genes=5000, 
    batch_key=snakemake.params.sample_key)

# Define mitochondria and ribosome genes to remove
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Make a copy of the AnnData atlas that only contains variable genes
filtered_adata = adata[:, (adata.var['highly_variable']) & ~(adata.var['mt']) & ~(adata.var['rb'])].copy()

# Setup SCVI on the data layer
scvi.model.SCVI.setup_anndata(
    filtered_adata, layer="data", batch_key=snakemake.params.sample_key)

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
    accelerator='cpu',  
    early_stopping=True,
    early_stopping_patience=20
)

# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(snakemake.output.model_history, index=False)

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
adata.write_h5ad(snakemake.output.merged_rna_anndata, compression='gzip')

model.save(snakemake.params.model, overwrite=True)
