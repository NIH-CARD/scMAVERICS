import snapatac2 as snap
import anndata as ad
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import episcanpy as epi
import decoupler as dc

# Load in the celltype atac object
cell_type_atac = sc.read_h5ad(snakemake.input.atac_anndata)

# Read in parameters
cell_type = snakemake.params.cell_type
disease_name = snakemake.params.disease
control_name = snakemake.params.control
disease_param = snakemake.params.disease_param

# Get the list of cells enriched for disease state
cell_type_control = cell_type_atac.obs[disease_param] == control_name
cell_type_disease = cell_type_atac.obs[disease_param] == disease_name

# Run differential expression test on all of the bins between control and given disease state
cell_type_diff_bins = snap.tl.diff_test(
    cell_type_atac,
    cell_group1 = cell_type_control,
    cell_group2 = cell_type_disease,
    covariates=['Age']
    direction = 'both',
    min_log_fc = 0
)

# Change the results to a DataFrame
cell_type_diff_df = pd.DataFrame(cell_type_diff_bins)
# Define the column names, reindex
cell_type_diff_df.columns = ['feature name','log2(fold_change)','p-value','adjusted p-value']
cell_type_diff_df.index = cell_type_diff_df['feature name']

# Define the chromosome as well as starts and stops
cell_type_diff_df['chr'] = [x.split(':')[0] for x in cell_type_diff_df['feature name']]
cell_type_diff_df['start'] = [int(x.split(':')[1].split('-')[0]) for x in cell_type_diff_df['feature name']]
cell_type_diff_df['end'] = [int(x.split(':')[1].split('-')[1]) for x in cell_type_diff_df['feature name']]

# Make sure the adjusted p-values are float types and log adjusted
cell_type_diff_df['adjusted p-value'] = cell_type_diff_df['adjusted p-value'].astype(float)
cell_type_diff_df['-log10(p-value)'] = -np.log10(cell_type_diff_df['adjusted p-value'])

# File save location
file_save = snakemake.output.output_DAR_data
cell_type_diff_df.to_csv(file_save)

# File save location
image_save = snakemake.output.output_figure

# Plot the volcano plot
dc.plot_volcano_df(
    cell_type_diff_df,
    x = 'log2(fold_change)',
    y = 'adjusted p-value',
    top = 20,
    lFCs_thr = .5,
    sign_thr = 1e-2,
    figsize = (4, 4),
    dpi = 600,
    return_fig = False,
    save = image_save
)
