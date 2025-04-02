import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import os
import anndata as ad
import scanpy as sc
from pycisTopic.cistopic_class import *

# Merge cistopic objects
cistopic_obj = merge([pickle.load(open(cistopic_path, "rb")) for cistopic_path in snakemake.input.cistopic_objects])

# Test if sample can be exported
pickle.dump(
    cistopic_obj,
    open(snakemake.output.merged_cistopic_object, "wb")
)

print('CisTopic object merged and saved')

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)
cell_data = rna.obs
cell_data['barcode'] = [x.split('_')[0] for x in cell_data.index]
# Add the sample_id variable
cell_data['sample_id'] = cell_data[snakemake.params.sample_key]

# Test if sample can be exported
pickle.dump(
    cistopic_obj,
    open(snakemake.output.merged_cistopic_object, "wb")
)

# Create AnnData object
adata = ad.AnnData(cistopic_obj.fragment_matrix.T)
adata.obs.index = cistopic_obj.cell_data['atlas_identifier'].to_list()
adata.var.index = cistopic_obj.region_names

# Add chromosome parameters
adata.var['chromosome'] = [x.split(':')[0] for x in adata.var.index]
adata.var['start'] = [x.split(':')[1].split('-')[0] for x in adata.var.index]
adata.var['end'] = [x.split(':')[1].split('-')[1] for x in adata.var.index]
adata.var['peak length'] = [int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) for x in adata.var.index]

transfer_params = [
    'atlas_identifier', 
    'Unique_nr_frag', 
    'cisTopic_nr_frag',
    'sample_id',
    'barcode',
    'cell_type', 
    'Primary Diagnosis']

# Create DataFrame of cisTopic sample parameters
cistopic_frag_data = cistopic_obj.cell_data[transfer_params].reset_index()
cistopic_frag_data.index = cistopic_frag_data['atlas_identifier']

# Add sample, barcode, cell type, and Primary diagnosis, number of fragments, and unique fragments to the anndata object
for param in transfer_params[1:]:
    barcode2param = cistopic_frag_data[param].to_dict()
    adata.obs[param] = [barcode2param[x] for x in adata.obs.index]

# Write out sample
adata.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')