import anndata as ad
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import muon as mu
import snapatac2 as snap

# Read in AnnData atlas object
mdata = mu.read(snakemake.input.multiome_object)

# Select for the most variable genes
sc.pp.highly_variable_genes(
    mdata.mod["rna"], 
    n_top_genes=snakemake.params.hvg,
    subset=False,
    flavor = 'seurat_v3')

rna_hvg = mdata["rna"][:, mdata["rna"].var['highly_variable']].copy()

snap.pp.select_features(mdata["atac"], n_features=snakemake.params.hvp, inplace=True)
atac_hvp = mdata["atac"][:, mdata["atac"].var["selected"].values].copy()

# Save the anndata object
mdata = mu.MuData({"rna": rna_hvg, "atac": atac_hvp})

mdata.write(snakemake.output.multiome_object)