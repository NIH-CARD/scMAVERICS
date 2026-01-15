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
    pl_fragment.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    # Get list of sample and cell type specific barcodes
    cell_type_barcodes = {cell_type: cell_df[(cell_df['celltype']== cell_type) & (cell_df[sample_value] == sample)]['cell_barcode'].to_list() for cell_type in snakemake.params.cell_types}
    # Filter on the cell type barcodes
    cell_fragment = {cell_type: pl_fragment.filter(pl_fragment['name'].is_in(cell_type_barcodes[cell_type])) for cell_type in snakemake.params.cell_types}
    # Add the filtered barcodes to the fragments
    print(f'Writing sample {sample}')
    for j, cell_type in enumerate(snakemake.params.cell_types):
        with open(snakemake.output.pseudo_fragment_files[j], mode='a') as f:
            cell_fragment[cell_type].write_csv(f, include_header=False, separator='\t')
            f.close()
    print(f'Sample {sample} has been added')

# Uncomment once singularity image is updated
#bigwig_path = f'/data/CARD_singlecell/PFC_atlas/data/celltypes/{cell_type}_test.bw'
#combined_bed.to_bigwig(bigwig_path, rpm=True)
