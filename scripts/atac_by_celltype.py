import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pycisTopic.cistopic_class import *
import scanpy as sc
import anndata as ad
import pyranges as pr

sample_key = snakemake.params.sample_key
seq_batch_key = snakemake.params.seq_batch_key
disease_param = snakemake.params.disease_param

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)
cell_data = rna.obs
print(rna.n_obs)

# Get cell barcodes of the given cell type
cell_data = cell_data[cell_data['cell_type'] == snakemake.params.cell_type]

cell_data['barcode'] = [x.split('_')[0] for x in cell_data.index]
cell_data['sample_id'] = cell_data[sample_key]
cell_sample_batch = cell_data[['sample_id', seq_batch_key]].drop_duplicates()

# Get fragment files and corresponding samples
fragment_files = snakemake.input.fragment_files
samples = snakemake.params.samples
fragment_dict = dict(zip(samples, fragment_files))

cell_type_samples = [x for x in cell_sample_batch['sample_id'].to_list()]
cell_type_fragments = {key: fragment_dict[key] for key in cell_type_samples}

adatas = []
print(f'Iterating through {len(cell_type_samples)} samples')
for sample, fragment_file in cell_type_fragments.items():
    cistopic_obj = create_cistopic_object_from_fragments(path_to_fragments=fragment_file,
                                               path_to_regions=snakemake.input.cell_bedfile,
                                               valid_bc = cell_data[cell_data['sample_id'] == sample]['barcode'].to_list(),
                                               n_cpu=1,
                                               project=sample
                                               )

    cistopic_obj.cell_data['atlas_identifier'] = [cistopic_obj.cell_data['barcode'][x] + '_' + cistopic_obj.cell_data['sample_id'][x] for x in range(len(cistopic_obj.cell_data))]
    adata = ad.AnnData(cistopic_obj.fragment_matrix.T)
    adata.obs = pd.merge(
        left=cistopic_obj.cell_data,
        right=cell_data[cell_data[sample_key] == sample][[disease_param, 'Age', seq_batch_key, 'PMI', 'atlas_identifier']],
        left_on='atlas_identifier',
        right_on='atlas_identifier')
    adata.var.index = cistopic_obj.region_names
    adata.var['chromosome'] = [x.split(':')[0] for x in adata.var.index]
    adata.var['start'] = [x.split(':')[1].split('-')[0] for x in adata.var.index]
    adata.var['end'] = [x.split(':')[1].split('-')[1] for x in adata.var.index]
    adata.var['peak length'] = [int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) for x in adata.var.index]
    adatas.append(adata)
    print(f"Done with sample {sample}")

adata = ad.concat(join='outer', adatas=adatas)

adata.var['chromosome'] = [x.split(':')[0] for x in adata.var.index]
adata.var['start'] = [x.split(':')[1].split('-')[0] for x in adata.var.index]
adata.var['end'] = [x.split(':')[1].split('-')[1] for x in adata.var.index]
adata.var['peak length'] = [int(x.split(':')[1].split('-')[1]) - int(x.split(':')[1].split('-')[0]) for x in adata.var.index]

adata.obs['atlas_identifier'] = adata.obs['atlas_identifier'].astype('category')
adata.obs_names = adata.obs['atlas_identifier'].astype(str).to_list()

# Convert the cisTopic parameters that were saved as objects?
adata.obs['cisTopic_nr_frag'] = adata.obs['cisTopic_nr_frag'].astype(int)
adata.obs['cisTopic_nr_acc'] = adata.obs['cisTopic_nr_acc'].astype(int)
adata.obs['atlas_identifier'] = adata.obs['atlas_identifier'].astype(str)

# Add required covariate values
for covariate in snakemake.params.covariates:
    cov_dict = cell_data[covariate].to_dict()
    adata.obs[covariate] = [cov_dict[x] for x in adata.obs_names.to_list()]

# Remove unnecessary outputs
del adata.obs['cisTopic_log_nr_frag']
del adata.obs['cisTopic_log_nr_acc']

annotated_bed = pd.read_csv(snakemake.input.cell_annotated_bedfile, delimiter='\t').iloc[:, 1:]
annotated_bed.index = annotated_bed['Chr'].astype(str) + ':' + (annotated_bed['Start']-1).astype(str) + '-' + annotated_bed['End'].astype(str)

adata.var = pd.merge(
    left=adata.var,
    right=annotated_bed[['Annotation', 'Detailed Annotation', 'Distance to TSS', 'Gene Name', 'Gene Type']],
    left_index=True,
    right_index=True)

adata.var['Annotation'] = adata.var['Annotation'].astype('category')
adata.var['Detailed Annotation'] = adata.var['Detailed Annotation'].astype('category')
adata.var['Gene Name  '] = adata.var['Gene Name'].astype('category')
adata.var['Gene Type '] = adata.var['Gene Type'].astype('category')
adata.var['start'] = adata.var['start'].astype(int)
adata.var['end'] = adata.var['end'].astype(int)
adata.var['total counts'] = np.array(adata.X.sum(axis=0))[0]
adata.var['Simple Annotation'] = [str(x).split(' ')[0] for x in adata.var['Annotation']]

adata.write_h5ad(snakemake.output.celltype_atac, compression='gzip')
