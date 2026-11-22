import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import decoupler as dc
import muon as mu
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Open multiome object
mdata = mu.read_h5mu(snakemake.input.mdata)

# Slice out and save region data
region = snakemake.params.region
mdata = mdata[mdata.obs[snakemake.params.region_param] == region].copy()
mdata.write_h5mu(snakemake.output.region_mdata, compression = 'gzip')

# Use only the rna modality
adata = mdata.mod['rna']

# Get pseudo-bulk profile, summing all raw counts by cell type and sample
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=snakemake.params.sample_key,
    groups_col=snakemake.params.separating_cluster,
    layer='counts',
    mode='sum'
)

dc.pp.filter_samples(pdata, min_cells=snakemake.params.min_cells)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

dc.pp.filter_by_expr(
    adata=pdata,
    group="celltype",
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
)
dc.pp.filter_by_prop(
    adata=pdata,
    min_prop=0.1,
    min_smpls=2,
)

pdata.write_h5ad(snakemake.output.pseudo_rna, compression='gzip')

work_dir+'atlas/STL_pseudobulked_rna.h5ad'
work_dir+'atlas/STL_pseudobulked_chromvar.h5ad'
work_dir+'atlas/STL_pseudobulked_atac.h5ad'