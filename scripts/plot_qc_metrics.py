import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np

# Keep consistent font sizes

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) 

# Convert to metric
cm = 1/2.54

# Set figure context
sns.set_context('paper')

# Define sample key
sample_key = snakemake.params.sample_key

# Open the AnnData object
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata) 

# Check for plots directory and create if not there
os.makedirs('figures/plots', exist_ok=True)

# Create dummy QC AnnData obj to concat doublet-calculated samples
doublet_adata = adata[adata.obs['total_counts'] < 0]

for sample in adata.obs[sample_key].drop_duplicates().to_list():

    # Make plot directory
    try:
        os.mkdir(f'plots/{sample}')
    except FileExistsError:
        print('Already there')


    """Plot percent mitochondria"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs[sample_key] == sample], ['pct_counts_mt'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('percent')
    ax[0].set_xticks('')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [snakemake.params.mito_percent_thresh, snakemake.params.mito_percent_thresh], '--r')
    ax[0].set_title('Percent mitochondria per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[adata.obs[sample_key] == sample].obs['pct_counts_mt'], 
        bins=int(np.sqrt(adata[adata.obs[sample_key] == sample].n_obs))
        )
    ax[1].set_yscale("log")
    ax[1].set_xlabel('percent')
    ax[1].set_ylabel('number of cells')
    ax[1].plot([snakemake.params.mito_percent_thresh, snakemake.params.mito_percent_thresh], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_title('Percent mitochondria per cell')
    plt.savefig(f'figures/plots/{sample}/mito_pct.png', dpi=300)


    """Plot percent ribosome"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs[sample_key] == sample], ['pct_counts_rb'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('percent')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [snakemake.params.ribo_percent_thresh, snakemake.params.ribo_percent_thresh], '--r')
    ax[0].set_title('Percent ribosome genes per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[adata.obs[sample_key] == sample].obs['pct_counts_rb'], 
        bins=int(np.sqrt(adata[adata.obs[sample_key] == sample].n_obs))
        )
    ax[1].plot([snakemake.params.ribo_percent_thresh, snakemake.params.ribo_percent_thresh], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_yscale("log")
    ax[1].set_xlabel('percent')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Percent ribosome genes per cell')
    plt.savefig(f'figures/plots/{sample}/ribo_pct.png', dpi=300)


    """Plot number of genes per cell"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs[sample_key] == sample], ['n_genes_by_counts'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('total counts')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [snakemake.params.min_genes_per_cell, snakemake.params.min_genes_per_cell], '--r')
    ax[0].set_title('Number of genes per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[adata.obs[sample_key] == sample].obs['n_genes_by_counts'], 
        bins=int(np.sqrt(adata[adata.obs[sample_key] == sample].n_obs))
        )
    ax[1].plot([snakemake.params.min_genes_per_cell, snakemake.params.min_genes_per_cell], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_xlabel('total counts')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Number of genes per cell')
    plt.savefig(f'figures/plots/{sample}/num_genes_per_cell.png', dpi=300)


    """Plot the scrublet values"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Run scrublet on filtered samples 
    filtered_adata = adata[
        (adata.obs[sample_key] == sample) & 
        (adata.obs['pct_counts_mt'] < snakemake.params.mito_percent_thresh) &
        (adata.obs['n_genes_by_counts'] > snakemake.params.min_genes_per_cell) &
        (adata.obs['pct_counts_rb'] < snakemake.params.ribo_percent_thresh)].copy()
    # Rerun scrublet
    sc.pp.scrublet(filtered_adata, expected_doublet_rate=(filtered_adata.n_obs / 1000) * 0.008, threshold=0.15, n_prin_comps=10)
    # Save the scrublet values to a new AnnData object
    doublet_adata = ad.concat([doublet_adata, filtered_adata], join='outer')

    # Violin plot in the first panel
    sc.pl.violin(filtered_adata, ['doublet_score'], jitter=0.5, ax=ax[0], show=False)
    ax[0].plot([-.5, .5], [snakemake.params.doublet_thresh, snakemake.params.doublet_thresh], '--r')
    ax[0].set_ylabel('droplet score')
    ax[0].set_xlim(-.75, .75)
    ax[0].set_title('Doublet score per cell')

    # Histogram
    y, x, _ = ax[1].hist(
        filtered_adata, 
        bins=int(filtered_adata.n_obs))
    ax[1].plot([snakemake.params.doublet_thresh, snakemake.params.doublet_thresh], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_xlabel('droplet score')
    ax[1].set_ylabel('number of droplets')
    ax[1].set_title('Doublet score per cell')
    plt.savefig(f'figures/plots/{sample}/scrublet_score_per_cell.png', dpi=300)

    # Plot number of genes by number of counts per cell
    sc.pl.scatter(
        adata[adata.obs[sample_key] == sample], 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_mt",
        show=False
        )
    plt.savefig(f'figures/plots/{sample}/num_gene_counts_total.png', dpi=300)


# Plot summary mitochondria, ribosome, and scrublet scores

# Mitochondria QC
y, x, _ = sns.hist(
    adata.obs.obs['pct_counts_mt'], 
    bins=int(np.sqrt(adata.n_obs))
    )
plt.yscale("log")
plt.xlabel('percent')
plt.ylabel('number of cells')
plt.plot([snakemake.params.mito_percent_thresh, snakemake.params.mito_percent_thresh], [1, y.max()], '--r')
plt.ylim(0, y.max())
plt.title('Percent mitochondria per cell')
plt.savefig(snakemake.output.mito_figure, dpi=300)

# Ribosome QC
y, x, _ = sns.hist(
        adata.obs['pct_counts_rb'], 
        bins=int(np.sqrt(adata.n_obs))
        )
plt.plot([snakemake.params.ribo_percent_thresh, snakemake.params.ribo_percent_thresh], [1, y.max()], '--r')
plt.ylim(0, y.max())
plt.yscale("log")
plt.xlabel('percent')
plt.ylabel('number of cells')
plt.title('Percent ribosome genes per cell')
plt.savefig(snakemake.output.ribo_figure, dpi=300)

# Number of genes
y, x, _ = sns.hist(
        adata.obs['n_genes_by_counts'], 
        bins=int(np.sqrt(adata.n_obs))
        )
plt.plot([snakemake.params.min_genes_per_cell, snakemake.params.min_genes_per_cell], [1, y.max()], '--r')
plt.ylim(0, y.max())
plt.xlabel('total counts')
plt.ylabel('number of cells')
plt.title('Number of genes per cell')
plt.savefig(snakemake.output.gene_counts_figure, dpi=300)

# Doublet QC
y, x, _ = sns.hist(
        doublet_adata, 
        bins=int(doublet_adata.n_obs))
plt.plot([snakemake.params.doublet_thresh, snakemake.params.doublet_thresh], [1, y.max()], '--r')
plt.ylim(0, y.max())
plt.xlabel('droplet score')
plt.ylabel('number of droplets')
plt.title('Doublet score per cell')
plt.savefig(snakemake.output.doublet_figure, dpi=300)

# Genes by total counts
sc.pl.scatter(
    adata, 
    "total_counts", 
    "n_genes_by_counts", 
    color="pct_counts_mt",
    show=False
    )
plt.savefig(snakemake.output.genes_by_counts, dpi=300)