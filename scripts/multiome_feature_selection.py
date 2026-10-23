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

# Read in ATAC dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Select for the most variable genes
sc.pp.highly_variable_genes(
    rna, 
    n_top_genes=snakemake.params.hvg,
    subset=False,
    flavor = 'seurat_v3')
# Make a copy that only contains the highly variable genes
rna_hvg = rna[:, rna.var['highly_variable']].copy()

# Select for the most variable peaks
snap.pp.select_features(atac, n_features=snakemake.params.hvp, inplace=True)
# Make a copy that only contains the highly variable peaks
atac_hvp = atac[:, atac.var["selected"].values].copy()

# Verify the RNA layer is using counts (this is important for MultiVI)
if "counts" in rna_hvg.layers:
    rna_hvg.X = rna_hvg.layers["counts"].copy()

# Save modality specific barcodes
rna_barcodes = rna_hvg.obs_names
atac_barcodes = atac_hvp.obs_names

# Find the paired, RNA-, and ATAC-only barcodes
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
    X=sp.csr_matrix((len(rna_only_barcodes), atac_hvp.n_vars), dtype=atac_hvp.X.dtype),
    obs=pd.DataFrame(index=rna_only_barcodes),
    var=atac_hvp.var.copy()
)

# Create RNA and ATAC objects that have the same index, missing values
# will be filled the 'rna/atac_missing' objects
unified_rna = ad.concat([rna_hvg, rna_missing], merge="same")
unified_atac = ad.concat([atac_missing, atac_hvp], merge="same")

# Ensure only the nuclei with paired barcodes are included
unified_rna = unified_rna[all_barcodes, :]
unified_atac = unified_atac[all_barcodes, :]

# Create the object
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

# Create sample key for batch-corrected MultiVI
mdata.obs[snakemake.params.sample_key] = [(x.split('-1_')[0]) for x in mdata.obs_names]
# Save object
mdata.write_h5mu(snakemake.output.multiome_object)
