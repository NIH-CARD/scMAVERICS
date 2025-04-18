import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pycisTopic.cistopic_class import *
import scanpy as sc
import pickle 
import os

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)
cell_data = rna.obs

# Add the sample_id variable
cell_data['barcode'] = [x.split('_')[0] for x in cell_data.index]
cell_data['sample_id'] = cell_data[snakemake.params.sample_key]
sample_batch = cell_data[[snakemake.params.sample_key, snakemake.params.seq_batch_key]].drop_duplicates()

# Path to regions
path_to_regions = snakemake.input.consensus_bed
# Create cistopic object
cistopic_obj = create_cistopic_object_from_fragments(path_to_fragments=snakemake.input.fragment_file,
                                               path_to_regions=path_to_regions,
                                               valid_bc = cell_data[cell_data['Sample'] == snakemake.params.sample]['barcode'].to_list(),
                                               n_cpu=32,
                                               project=snakemake.params.sample)

# Assign metadata
cistopic_obj.cell_data['atlas_identifier'] = [cistopic_obj.cell_data['barcode'][i] + '_' + cistopic_obj.cell_data['sample_id'][i] for i in range(len(cistopic_obj.cell_data))]

barcode2celltype = cell_data['cell_type'].to_dict()
barcode2disease = cell_data[snakemake.params.disease_param].to_dict()

cistopic_obj.cell_data['cell_type'] = [barcode2celltype[x] for x in cistopic_obj.cell_data['atlas_identifier']]
cistopic_obj.cell_data[snakemake.params.disease_param] = [barcode2disease[x] for x in cistopic_obj.cell_data['atlas_identifier']]

# Export sample
pickle.dump(
    cistopic_obj,
    open(snakemake.output.cistopic_object, "wb")
)

# Create AnnData object
adata = ad.AnnData(cistopic_obj.fragment_matrix.T)
adata.obs.index = cistopic_obj.cell_data['atlas_identifier'].to_list()
adata.var.index = cistopic_obj.region_names

# Add chromosome variable data
adata.var['chromosome'] = [x.split(':')[0] for x in adata.var.index]
adata.var['start'] = [x.split(':')[1].split('-')[0] for x in adata.var.index]
adata.var['end'] = [x.split(':')[1].split('-')[1] for x in adata.var.index]
adata.var['peak length'] = [int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) for x in adata.var.index]

# Variables to transfer
transfer_params = [
    'atlas_identifier', 
    'Unique_nr_frag', 
    'cisTopic_nr_frag',
    'sample_id',
    'barcode',
    'cell_type', 
    snakemake.params.disease_param]

# Convert variable metadata to DataFrame to be transfered 
cistopic_frag_data = cistopic_obj.cell_data[transfer_params].reset_index()
cistopic_frag_data.index = cistopic_frag_data['atlas_identifier']

# Add sample, barcode, cell type, and Primary diagnosis, number of fragments, and unique fragments to the anndata object
for param in transfer_params[1:]:
    barcode2param = cistopic_frag_data[param].to_dict()
    adata.obs[param] = [barcode2param[x] for x in adata.obs.index]

# Write out peak AnnData object
adata.write_h5ad(snakemake.output.cistopic_adata, compression='gzip')