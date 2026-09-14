import numpy as np
import pandas as pd
import scanpy as sc
import polars as pl
import pyranges as pr

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Read in atac observation data
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)
# Filter so only barcodes in both samples are used for creating fragments
rna = rna[rna.obs_names.intersection(atac.obs_names)]

# Port cell data from final RNA atlas to cisTopic pseudobulked
cell_df = rna.obs

if 'cell_barcode' not in cell_df.columns:
    cell_df['cell_barcode'] = [x.split('-1_')[0] + '-1' for x in cell_df.index]

# Metadata specific column names
sample_value = snakemake.params.sample_param_name

# Get sample list
samples = snakemake.params.samples

# For sample in samples
for i, sample in enumerate(samples):
    # Define fragment file
    bed_location = snakemake.input.fragment_file[i]
    # Load fragment with polars
    print(f'Loading sample {sample} fragments')
    pl_fragment = pl.read_csv(bed_location, separator='\t', comment_prefix='#', n_threads=8)

    # To make sure the same columns are being read in
    num_blank = len(pl_fragment.columns) - 5
    pl_fragment.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score'] + ['.']*num_blank
    pl_fragment.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score'] + ['.']*num_blank
    pl_fragment = pl_fragment[['chrom', 'chromStart', 'chromEnd', 'name', 'score']]
    
    # Get list of sample and cell type specific barcodes
    cell_type_barcodes = cell_df[(cell_df['celltype']== snakemake.params.cell_type) & (cell_df[sample_value] == sample)]['cell_barcode'].to_list()
    # Filter on the cell type barcodes
    cell_fragment = pl_fragment.filter(pl_fragment['name'].is_in(cell_type_barcodes))
    # Add the filtered barcodes to the fragments
    print(f'Writing sample {sample}')
    with open(snakemake.output.pseudo_fragment_files, mode='a') as f:
        cell_fragment.write_csv(f, include_header=False, separator='\t')
        f.close()
    print(f'Sample {sample} has been added')
