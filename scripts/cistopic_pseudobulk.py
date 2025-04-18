import numpy as np
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
import scanpy as sc
import os

# Read in rna observation data
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Port cell data from final RNA atlas to cisTopic pseudobulked
cell_data = rna.obs
# Add the sample_id and barcode variables to match required cisTopic input
cell_data['barcode'] = [x.split('_')[0] for x in cell_data.index]
cell_data['sample_id'] = cell_data[snakemake.params.sample_param_name]

# Make sure list of samples is interpreted as strings
samples = [str(x) for x in snakemake.params.samples]

# Subset metadata to just samples in dataset 
cell_data = cell_data[cell_data['sample_id'].isin(samples)]

# Load chromosome sizes
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

# Input the fragment files with the same input
fragments_dict = dict(zip(samples, snakemake.input.fragment_file))


# Create folds for pseudobulked samples
os.makedirs(snakemake.params.bed_file_locs, exist_ok=True)
os.makedirs(snakemake.params.bigwig_file_locs, exist_ok=True)

# Export normalized pseudobulk bed and bigwig files
bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = snakemake.params.pseudobulk_param,
    chromsizes = chromsizes,
    bed_path = snakemake.params.bed_file_locs,
    bigwig_path = snakemake.params.bigwig_file_locs,
    path_to_fragments = fragments_dict,
    n_cpu = snakemake.threads,
    normalize_bigwig = False,
    temp_dir = snakemake.resources.tmpdir
    )
