import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import os
import anndata as ad
import scanpy as sc


# Load samples matched to 
sample_loc = dict(zip(snakemake.params.samples, snakemake.input.atac_anndata))

adatas = {key: sc.read_h5ad(sample_loc[key]) for key in sample_loc.keys()}

adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas=adatas
    )

# Write out sample
adata.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')