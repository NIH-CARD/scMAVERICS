import anndata as ad
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import scipy

# Read in anndata objects using backed mode
adatas = [snap.read(f, backed="r") for f in snakemake.input.adatas]
print('Done reading adatas')

# 
anndataset = snap.concat(
    adatas, 
    keys=snakemake.params.samples, 
    label=snakemake.params.sample_key,
    file=snakemake.output.merged_atac_anndata
)
print('Dataset merged')

# Populate merged object with calculated metadata
anndataset.obs_names = [
    f"{barcode}_{sample}" 
    for barcode, sample in zip(anndataset.obs_names, anndataset.obs[snakemake.params.sample_key])
]
anndataset.obs['n_fragment'] = anndataset.adatas.obs['n_fragment']
anndataset.obs['tsse'] = anndataset.adatas.obs['tsse']

# Close input files
for adata in adatas:
    adata.close()
print('Anndata objects closed')

# Close merged output file
anndataset.close()
print('Merged anndata object closed')