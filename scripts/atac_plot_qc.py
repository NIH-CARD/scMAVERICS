import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

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
adata = sc.read_h5ad(snakemake.input.merged_atac_anndata) 

# Check for plots directory and create if not there
os.makedirs('plots', exist_ok=True)

# Iterate through the samples
for sample in adata.obs[sample_key].to_list():

    # Make plot directory
    try:
        os.mkdir(f'plots/{sample}')
    except FileExistsError:
        print('Already there')

    """Plot percent mitochondria"""

    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    sc.pl.violin(
        adata[adata.obs[sample_key] == sample], 
        ['n_genes_by_counts'], 
        jitter=0.5, 
        ax=ax[0],
         show=False)
    ax[0].set_ylabel('number of open features')
    ax[0].set_xticks('')
    ax[0].set_xlim(-.5, .5)
    ax[0].plot([-.5, .5], [2000, 2000], '--r')
    ax[0].set_title('Number of open features per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[atac.obs[sample_key] == sample].obs['n_genes_by_counts'], 
        bins=int(np.sqrt(adata[adata.obs[sample_key] == sample].n_obs))
        )
    ax[1].set_xlabel('number of open features')
    ax[1].set_ylabel('number of cells')
    ax[1].plot([2000, 2000], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_title('Number of open features per cell')
    plt.savefig(f'plots/{sample}/atac_open_features_per_cell.png', dpi=300)


    """"""
    # Define the sample
    sample = adata.obs[sample_key].iloc[0]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    shared_fragments = np.array(adata.var['n_cells_by_counts'])
    sns.violinplot(np.log10(shared_fragments), ax=ax[0])
    sns.stripplot(
        np.log10(shared_fragments), 
        color='black', 
        alpha=.25, 
        jitter=.5, 
        size=1, 
        ax=ax[0]
        )
    ax[0].set_ylabel('number of shared open features')
    ax[0].set_xticks('')
    ax[0].set_xlim(-.5, .5)
    cell_by_count_thresh = np.log10(10)
    ax[0].plot([-.5, .5], [cell_by_count_thresh, cell_by_count_thresh], '--r')
    ax[0].set_title('Number of shared open features per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        np.log10(shared_fragments), 
        bins=int(np.sqrt(adata[adata.obs[sample_key] == sample].n_obs))
        )
    ax[1].set_xlabel('number of shared open features')
    ax[1].set_ylabel('number of cells')
    ax[1].plot([cell_by_count_thresh, cell_by_count_thresh], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_title('Number of shared open features per cell')
    plt.savefig(f'plots/{sample}/atac_shared_features_between_cell.png', dpi=300)
