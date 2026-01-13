import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import gseapy as gp

# Load in cell type AnnData object
adata = sc.read_h5ad(snakemake.input.adata_path)

adata = adata[adata.obs[snakemake.params.cell_param] == snakemake.params.cell_type]

# Filter based on condition to test
diseaes_adata = adata[adata.obs[snakemake.params.disease_param].isin([snakemake.params.control, snakemake.params.disease])]

# Add one-hot column for disease comparison, sort the AnnData object
diseaes_adata.obs['disease'] = pd.Categorical(
    diseaes_adata.obs[snakemake.params.disease_param], 
    categories=[snakemake.params.disease, snakemake.params.control], 
    ordered=True)
indices = diseaes_adata.obs.sort_values(['disease']).index
diseaes_adata = diseaes_adata[indices,:].copy()

# Import same ontology dataset as used in GREAT, reformat for GSEA
greatpy_go_set_df = pd.read_csv(snakemake.input.ontologies, delimiter=';')
greatpy_go_set_df['key'] = greatpy_go_set_df['name'] + ' (' + greatpy_go_set_df['id'] + ')'
go_mf = greatpy_go_set_df[['key', 'symbol']].groupby('key')['symbol'].apply(list).to_dict()

# Calculate GSEA enrichment
res = gp.gsea(
    data=diseaes_adata.to_df().T, # row -> genes, column-> samples
    gene_sets = go_mf,
    classes=[snakemake.params.disease, snakemake.params.control],
    permutation_num=1000,
    permutation_type='phenotype',
    outdir=None,
    method='log2_ratio_of_classes', # signal_to_noise
    threads= snakemake.threads # Speed up enrichment analysis
    )

# Save file
res.res2d.to_csv(snakemake.output.cell_disease_GSEA, index=False)