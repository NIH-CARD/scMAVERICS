import snapatac2 as snap
import pandas as pd
import numpy as np
import scanpy as sc


# Read list of atac data locations
atac_anndata = snakemake.input.atac_anndata
print(len(atac_anndata))

# Read in snapATAC2 datasets into a list of anndata objects in read only
list_of_anndata = [(samples_to_keep[i], snap.read(atac_anndata[i])) for i in range(len(atac_anndata))]

anndataset = snap.AnnDataSet(
    adatas=list_of_anndata,
    filename=snakemake.output.merged_atac_anndata,
    backed=True
)

# Update identifiers
anndataset.obs_names = [bc + '_' + sa for bc, sa in zip(anndataset.obs_names, anndataset.obs['sample'])]

# Add metadata
metadata_dict = {
    'batch' : 'Use_batch',
    'sex': 'Sex',
    'age': 'Age',
    'pmi': 'PMI',
    'ethnicity': 'Ethnicity',
    'race': 'Race',
    'brain_bank': 'Brain_bank',
    'diagnosis': 'PrimaryDiagnosis'
}

for new_obs, metadata_obs in metadata_dict.items():
    anndataset.obs[new_obs] = [samples.loc[samples['Sample_ID'] == x][metadata_obs].iloc[0] for x in anndataset.obs['sample']]

# Select variable features
snap.pp.select_features(adataset, n_features=250000, n_jobs=60)

# Spectral MDS analysis
snap.tl.spectral(adataset)

# Batch correction
snap.pp.mnc_correct(adataset, batch="sample", key_added='X_spectral')

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
rna_annot = pd.read_csv('/data/CARD_singlecell/SN_atlas/data/rna_cell_annot.csv')
anndataset.obs['cell_type'] = rna_annot['cell_type'].to_list()

# Close dataset
anndataset.close()
