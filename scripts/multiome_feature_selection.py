import anndata as ad
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np
import muon as mu
import snapatac2 as snap

# Read in RNA dataset
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Select for the most variable genes
sc.pp.highly_variable_genes(
    rna, 
    n_top_genes=snakemake.params.hvg,
    subset=False,
    flavor = 'seurat_v3')

rna_hvg = rna[:, rna.var['highly_variable']].copy()

# Read in ATAC dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

snap.pp.select_features(mdata["atac"], n_features=snakemake.params.hvp, inplace=True)
atac_hvp = mdata["atac"][:, mdata["atac"].var["selected"].values].copy()

# Use raw RNA counts
if "counts" in rna.layers:
    rna_hvg.X = rna_hvg.layers["counts"].copy()

# Save modality specific barcodes
rna_barcodes = rna_hvg.obs_names
atac_barcodes = atac_hvp.obs_names

all_barcodes = rna_barcodes.union(atac_barcodes, sort=False)
paired_barcodes = rna_barcodes.intersection(atac_barcodes, sort=False)
rna_only_barcodes = rna_barcodes.difference(atac_barcodes, sort=False)
atac_only_barcodes = atac_barcodes.difference(rna_barcodes, sort=False)

# Keep original RNA object; only create missing rows for ATAC-only nuclei
rna_missing = ad.AnnData(
    X=sp.csr_matrix((len(atac_only_barcodes), rna_hvg.n_vars), dtype=rna_hvg.X.dtype),
    obs=pd.DataFrame(index=atac_only_barcodes),
    var=rna_hvg.var.copy()
)

# Keep original ATAC object; only create missing rows for RNA-only nuclei
atac_missing = ad.AnnData(
    X=sp.csr_matrix((len(rna_only_barcodes), atac_hvp.n_vars), dtype=atac.X.dtype),
    obs=pd.DataFrame(index=rna_only_barcodes),
    var=atac_hvp.var.copy()
)

unified_rna = ad.concat([rna_hvg, rna_missing], merge="same")
unified_atac = ad.concat([atac_missing, atac_hvp], merge="same")

unified_rna = unified_rna[all_barcodes, :]
unified_atac = unified_atac[all_barcodes, :]

mdata = mu.MuData({"rna": unified_rna, "atac": unified_atac})

# Metadata for plotting
paired_set = set(paired_barcodes)
rna_only_set = set(rna_only_barcodes)
atac_only_set = set(atac_only_barcodes)

mdata.obs["modality_group"] = pd.Categorical(
    [
        "paired" if x in paired_set else
        "rna_only" if x in rna_only_set else
        "atac_only"
        for x in mdata.obs_names
    ],
    categories=["paired", "rna_only", "atac_only"]
)

mdata.write(snakemake.output.hvg_multiome_anndata)