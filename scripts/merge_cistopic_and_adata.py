import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import os
import anndata as ad
import scanpy as sc

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)
cell_data = rna.obs
# Add the sample_id variable

samples = cell_data['SampleID'].to_list()

sample_loc = dict(zip(samples, snakemake.input.atac_anndata))

adatas = {key: sc.read_h5ad(sample_loc[key]) for key in samples}

adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas=adatas
    )

# Write out sample
adata.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')