import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import scipy.sparse as sp


# Read in RNA dataset
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Read in ATAC dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Use raw RNA counts
if "counts" in rna.layers:
    rna.X = rna.layers["counts"].copy()

# Save modality specific barcodes
rna_barcodes = rna.obs_names
atac_barcodes = atac.obs_names

all_barcodes = rna_barcodes.union(atac_barcodes, sort=False)
paired_barcodes = rna_barcodes.intersection(atac_barcodes, sort=False)
rna_only_barcodes = rna_barcodes.difference(atac_barcodes, sort=False)
atac_only_barcodes = atac_barcodes.difference(rna_barcodes, sort=False)

# Keep original RNA object; only create missing rows for ATAC-only nuclei
rna_missing = ad.AnnData(
    X=sp.csr_matrix((len(atac_only_barcodes), rna.n_vars), dtype=rna.X.dtype),
    obs=pd.DataFrame(index=atac_only_barcodes),
    var=rna.var.copy()
)

# Keep original ATAC object; only create missing rows for RNA-only nuclei
atac_missing = ad.AnnData(
    X=sp.csr_matrix((len(rna_only_barcodes), atac.n_vars), dtype=atac.X.dtype),
    obs=pd.DataFrame(index=rna_only_barcodes),
    var=atac.var.copy()
)

unified_rna = ad.concat([rna, rna_missing], merge="same")
unified_atac = ad.concat([atac_missing, atac], merge="same")

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
mdata.obs[snakemake.params.sample_key] = ['-1'.join(x.split('-1_')[1:]) for x in mdata.obs_names]

mdata.write(snakemake.output.multiome_object, compression ='gzip')
