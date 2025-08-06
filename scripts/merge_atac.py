import snapatac2 as snap
import pandas as pd
import numpy as np
import scanpy as sc

# Import the samples
samples = snakemake.params.samples
atac_anndata = snakemake.input.atac_anndata

# Import metadata
metadata_table = pd.read_csv(snakemake.input.metadata_table)
# Reset the index for the samples
metadata_table = metadata_table.set_index(snakemake.params.sample_key)

# Read in snapATAC2 datasets into a list of anndata objects in read only
list_of_anndata = [(samples[i], snap.read(atac_anndata[i])) for i in range(len(atac_anndata))]

# Create the AnnDataSet from the snapATAC2 datasets
anndataset = snap.AnnDataSet(
    adatas=list_of_anndata,
    # filename=atac_anndata
    filename=snakemake.output[0],
)

# Update identifiers
anndataset.obs_names = [bc + '_' + sa for bc, sa in zip(anndataset.obs_names, anndataset.obs[snakemake.params.sample_key])]

for value in metadata_table.columns.to_list():
    # Create a new dictionary for each sample-metadata value
    sample2value = metadata_table[value].to_dict()
    # Assign this value
    anndataset.obs[new_obs] = [sample2value[x] for x in anndataset.obs[snakemake.params.sample_key].to_list()]

# Select variable features
snap.pp.select_features(anndataset, n_features=250000, n_jobs=60)

# Spectral MDS analysis
snap.tl.spectral(anndataset)

# Batch correction
snap.pp.mnc_correct(anndataset, batch=snakemake.params.sample_key, key_added='X_spectral')

# Compute nearest neighbors from the corrected spectral MDS
snap.pp.knn(anndataset)

# Compute clusters
snap.tl.leiden(anndataset)

# Compute umap
snap.tl.umap(anndataset)

# Save values
pd.DataFrame(anndataset.obsm['X_umap']).to_csv(snakemake.output.umap_data)
pd.DataFrame(anndataset.var[['count', 'selected']]).to_csv(snakemake.output.var_data)

""" THIS AREA FOR INTEGRATING ANNOTATION WITH RNA DATA"""
# Save the dataframe
rna_annot = pd.read_csv(snakemake.input.cell_annotate)
anndataset.obs['cell_type'] = rna_annot['cell_type'].to_list()

# Close dataset
anndataset.close()
