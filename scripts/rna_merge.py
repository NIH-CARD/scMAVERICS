import anndata as ad
import pandas as pd
import scanpy as sc
import scipy

# Create dictionary of sample name to file, then filter with relevant diagnosis
sample_loc = dict(zip(snakemake.params.samples, snakemake.input.rna_anndata))
adatas = {key: sc.read_h5ad(sample_loc[key]) for key in snakemake.params.samples}

# Concatenate with names of sample
adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas=adatas
    )

# Produce a sparce counts layer
adata.layers['counts'] = scipy.sparse.csr_matrix(adata.layers['counts'].copy())
adata.layers['cpm'] = scipy.sparse.csr_matrix(adata.layers['cpm'].copy())
adata.layers['log-norm'] = scipy.sparse.csr_matrix(adata.layers['log-norm'].copy())

# Write out the unfiltered dataset
adata.write_h5ad(
    filename=snakemake.output.merged_rna_anndata, 
    compression='gzip'
    )