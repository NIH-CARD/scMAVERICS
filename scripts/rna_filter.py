import scanpy as sc

# Load in data
adata = sc.read_h5ad(snakemake.input.rna_anndata) # type: ignore

# Filter for quality control on general values
# The lambda function sets the observation row x to be filtered the variable in 
# quotes by the value input by the snakemake params value

# Threshold below a given mitochondria percent
adata = adata[adata.obs['pct_counts_mt'] < snakemake.params.mito_percent_thresh].copy()
# Threshold above a given genes per cell threshold
adata = adata[adata.obs['total_counts'] > snakemake.params.min_genes_per_cell].copy()
# Threshold below a given ribosome threshold
adata = adata[adata.obs['pct_counts_rb'] < snakemake.params.ribo_percent_thresh].copy()

# Scrublet requires raw counts
adata.X = adata.layers['counts']
# Compute doublet score based on filtered data
sc.pp.scrublet(adata, expected_doublet_rate=(adata.n_obs / 1000) * 0.008, threshold=0.15, n_prin_comps=10)
# Threshold below a given doublet score
adata = adata[adata.obs['doublet_score'] < snakemake.params.doublet_thresh].copy()

# Recompute log-normalized data from filtered dataset
if adata.n_obs != 0:
    # Normalize data
    sc.pp.normalize_total(adata)

    # Save the CPM data
    adata.layers['cpm']=adata.X.copy() 

    # Logarithmize the data
    sc.pp.log1p(adata, layer='cpm')

    # Save the normalized-log data
    adata.layers['log-norm']=adata.X.copy() 

# Write out filtered anndata object
adata.write(filename=snakemake.output.rna_anndata, compression='gzip') # type: ignore
