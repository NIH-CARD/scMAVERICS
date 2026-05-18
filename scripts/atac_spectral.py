import snapatac2 as snap
import scanpy as sc

# Read in merged dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Select the top N number of features 
snap.pp.select_features(atac, n_features=snakemake.params.num_features)

# Spectral multidimensional scaling
snap.tl.spectral(atac, num_threads=snakemake.threads)

# Perform mutual nearest neighbors
snap.pp.mnc_correct(
    atac, 
    batch=snakemake.params.sample_param, 
    key_added='X_spectral',
    n_jobs = snakemake.threads)

# Compute 
snap.pp.knn(
    atac, 
    use_rep='X_spectral'
    )

# Compute UMAP
snap.tl.umap(
    atac, 
    n_comps=2)

# Save file
atac.write_h5ad(snakemake.output.merged_atac_anndata, compression='gzip')