import numpy as np
import pandas as pd
import scanpy as sc
import polars as pl
import pyranges as pr

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Port cell data from final RNA atlas to cisTopic pseudobulked
cell_df = rna.obs

# Metadata specific column names
sample_key = snakemake.params.sample_param_name
disease_param = snakemake.params.disease_param

# Get snakemake params
samples = snakemake.params.samples
disease = snakemake.params.disease
cell_type = snakemake.params.cell_type

# Filter samples list
sample_list = cell_df[cell_df[disease_param] == disease][sample_key].drop_duplicates().to_list()

# For sample in samples
for i, sample in enumerate(samples):
    # Filter out samples read in, this should be replaced at some point with 
    # a smaller input fragment/sample lists
    if sample in sample_list:
        # Define fragment file
        bed_location = snakemake.input.fragment_file[i]
        # Load fragment with polars
        print(f'Loading sample {sample} fragments')
        pl_fragment = pl.read_csv(bed_location, separator='\t', comment_prefix='#', n_threads=8)
        pl_fragment.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
        # Get list of sample and cell type specific barcodes
        barcodes = cell_df[(cell_df['celltype'] == cell_type) & (cell_df[sample_key] == sample)]['cell_barcode'].to_list()

        # Filter on the cell type barcodes
        cell_fragment = pl_fragment.filter(pl_fragment['name'].is_in(barcodes))
        # Add the filtered barcodes to the fragments
        print(f'Writing sample {sample}')
    
        with open(snakemake.output.pseudo_fragment_files, mode='a') as f:
            cell_fragment[cell_type].write_csv(f, include_header=False, separator='\t')
            f.close()
