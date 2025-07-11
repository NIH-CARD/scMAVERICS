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

# Initiate polars dataframe
pseudo_fragments = pl.DataFrame()

# Filter for only the cell type data
cell_type_df = cell_df[cell_df[snakemake.params.pseudobulk_param] == snakemake.params.cell_type]

# Get sample list
samples = snakemake.params.samples
#batches = cell_type_df[[sample_value, batch_value]].drop_duplicates()[batch_value]

# For sample in samples
for i, sample in enumerate(samples):
    # Define fragment
    bed_location = snakemake.input.fragment_file[i]
    # Load fragment with polars
    print(f'Loading sample {sample} fragments')
    pl_fragment = pl.read_csv(bed_location, separator='\t', comment_prefix='#', n_threads=32)
    pl_fragment.columns = ['chr', 'start', 'end', 'barcode', 'count']
    # Get list of sample and cell type specific barcodes
    cell_type_barcodes = cell_type_df[(cell_type_df[sample_value] == sample)]['cell_barcode'].to_list()
    # Filter on the cell type barcodes
    pl_fragment = pl_fragment.filter(pl_fragment['barcode'].is_in(cell_type_barcodes))
    # Add the filtered barcodes to the fragments
    pseudo_fragments = pl.concat([pseudo_fragments, pl_fragment])
    print(f'Sample {sample} has been added')
    # Remove old dataframe
    del pl_fragment

# Filter only for chromosomes
filtered_fragments = pseudo_fragments.filter(pseudo_fragments['chr'].is_in(pr.data.chromsizes().df['Chromosome'].to_list()))
# Convert to pyranges
combined_bed = pr.PyRanges(
    chromosomes = filtered_fragments['chr'],
    starts=filtered_fragments['start'],
    ends=filtered_fragments['end']
    )

bed_path = f'/data/CARD_singlecell/SN_atlas/data/celltypes/{snakemake.params.cell_type}_bigwig.bed'
combined_bed.to_bed(bed_path, keep=False)

# Uncomment once singularity image is updated
bigwig_path = f'/data/CARD_singlecell/SN_atlas/data/celltypes/{snakemake.params.cell_type}_normalized_bigwig.bw'
combined_bed.to_bigwig(bigwig_path, rpm=True)
