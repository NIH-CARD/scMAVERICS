import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np

scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')

# Read in AnnData atlas object
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Produce a sparce counts layer
adata.layers['counts'] = scipy.sparse.csr_matrix(adata.layers['counts'].copy())
adata.layers['cpm'] = scipy.sparse.csr_matrix(adata.layers['cpm'].copy())
adata.layers['data'] = scipy.sparse.csr_matrix(adata.layers['data'].copy())

# Double check that no transcripts not found in cells are in the atlas
sc.pp.filter_genes(adata, min_cells=1)

# Select for the most variable genes
sc.pp.highly_variable_genes(
    adata, 
    layer='data',
    n_top_genes=5000, 
    batch_key="sample")

# Define mitochondria and ribosome genes to remove
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Make a copy of the AnnData atlas that only contains variable genes
filtered_adata = adata[:, (adata.var['highly_variable']) & ~(adata.var['mt']) & ~(adata.var['rb'])].copy()

# Setup SCVI on the data layer
scvi.model.SCVI.setup_anndata(
    filtered_adata, layer="data", batch_key="sample")

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
    accelerator='gpu',  
    early_stopping=True,
    early_stopping_patience=20
)

# Extract the elbo plot of the model and save the values
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(snakemake.output.model_history, index=False)

# 
adata.obs['atlas_identifier'] = adata.obs.index.to_list()
adata.obsm['X_scvi'] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=2, key_added='leiden_2')
sc.tl.leiden(adata, key_added='leiden')
sc.tl.leiden(adata, resolution=.5, key_added='leiden_05')

adata.write_h5ad(snakemake.output.merged_rna_anndata, compression='gzip')

model.save(snakemake.params.model, overwrite=True)
