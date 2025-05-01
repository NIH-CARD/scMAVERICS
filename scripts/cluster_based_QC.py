import sys
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

sc.pl.umap(
    adata,
    color='cell_type',
    legend_loc='on data',
    show=False
)
plt.savefig(snakemake.output.course_celltype)

# Plot and save the DataFrame to
sns.displot(
    adata.obs, 
    x="n_genes_by_counts",
    row="cell_type",
    binwidth=10, 
    hue='cell_type',
    height=.75, 
    aspect=4,
    kde=True,
    facet_kws=dict(sharey=False)).set(title='')
plt.xlim(0, 10000)
plt.savefig(snakemake.output.course_counts)

cell_type_min = {}
for cell_type in adata.obs['cell_type'].drop_duplicates():
    print(f'Working on cell type {cell_type}')

    cell_counts = adata[adata.obs['cell_type'] == cell_type].obs['n_genes_by_counts']
    # Return minimum values
    kde = scipy.stats.gaussian_kde(cell_counts)
    x_range = np.linspace(min(cell_counts), max(cell_counts), 1000)
    y_range = kde(x_range)

    # Find and assign the minimum values of the distribution
    min_index = scipy.signal.argrelextrema(y_range, np.less)
    minimum_x, minimum_y = x_range[min_index],  y_range[min_index]

    # Find and assign the maximum values of the distribution
    max_index = scipy.signal.argrelextrema(y_range, np.greater)
    maximum_x, maximum_y = x_range[max_index],  y_range[max_index]

    # Assign value if the local minimum meets a given threshold]
    thresh_min_x = minimum_x[minimum_x < 2000]
    # Only assign a value if there are values
    if thresh_min_x.shape[0] != 0:
        cell_type_min[cell_type] = thresh_min_x[0]
    else:
        cell_type_min[cell_type] = 0

bad_cells = sum([(adata.obs['cell_type'] == x) & (adata.obs['n_genes_by_counts'] < cell_type_min[x]) for x in cell_type_min.keys()])
bad_bools = [True if x > 0 else False for x in bad_cells]

adata.obs['bad_cells'] = bad_bools
adata = adata[~adata.obs['bad_cells']].copy()

# Save file
adata.write_h5ad(snakemake.output.merged_rna_anndata)