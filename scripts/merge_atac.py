import snapatac2 as snap
import pandas as pd
import numpy as np
import scanpy as sc

# Import the samples
samples = snakemake.params.samples
atac_anndata = snakemake.input.atac_anndata

# Import metadata
metadata = pd.read_csv(snakemake.input.metadata_table)
# Reset the index for the samples
metadata = metadata_table.set_index(snakefile.params.sample_key)

# Read in snapATAC2 datasets into a list of anndata objects in read only
list_of_anndata = [(samples[i], snap.read(atac_anndata[i])) for i in range(len(atac_anndata))]

# Create the AnnDataSet from the snapATAC2 datasets
anndataset = snap.AnnDataSet(
    adatas=list_of_anndata,
    filename=atac_anndata
)

# Update identifiers
anndataset.obs_names = [bc + '_' + sa for bc, sa in zip(anndataset.obs_names, anndataset.obs[snakemake.params.sample_key])]

for value in metadata_table.columns.to_list():
    # Create a new dictionary for each sample-metadata value
    sample2value = metadata_table[value].to_dict()
    # Assign this value
    anndataset.obs[new_obs] = [sample2value[x] for x in anndataset.obs[snakefile.params.sample_key].to_list()]

# Select variable features
snap.pp.select_features(adataset, n_features=250000, n_jobs=60)

# Spectral MDS analysis
snap.tl.spectral(adataset)

# Batch correction
snap.pp.mnc_correct(adataset, batch=snakemake.params.sample_key, key_added='X_spectral')

# Compute nearest neighbors from the corrected spectral MDS
snap.pp.knn(adataset)

# Compute clusters
snap.tl.leiden(adataset)

# Compute umap
snap.tl.umap(adataset)

# Save values
pd.DataFrame(adataset.obsm['X_umap']).to_csv(snakemake.output.umap_data)
pd.DataFrame(adataset.var[['count', 'selected']])to_csv(snakemake.output.var_data)

""" THIS AREA FOR INTEGRATING ANNOTATION WITH RNA DATA"""
# Save the dataframe
rna_annot = pd.read_csv(snakemake.input.cell_annotate)
anndataset.obs['cell_type'] = rna_annot['cell_type'].to_list()

# Close dataset
anndataset.close()
