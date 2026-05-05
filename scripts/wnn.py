import scanpy as sc
import muon as mu

# Read in multiome object
mdata = mu.read(snakemake.input.merged_multiome)

# Compute atac cosine-knn
sc.pp.neighbors(
    mdata.mod['atac'],
    n_neighbors=snakemake.params.num_neighbors,
    metric='cosine', 
    use_rep=snakemake.params.rna_rep
)
sc.tl.umap(mdata.mod['atac'], min_dist=0.3)

# Compute rna cosine-knn
sc.pp.neighbors(
    mdata.mod['rna'], 
    n_neighbors=snakemake.params.num_neighbors,
    metric='cosine', 
    use_rep=snakemake.params.atac_rep
)
sc.tl.umap(mdata.mod['atac'], min_dist=0.3)

# Compute merged modality cosine-knn
mu.pp.neighbors(
    mdata, 
    key_added='wnn', 
    metric='cosine', 
    low_memory=True)

# Compute merged modality UMAP
mu.tl.umap(
    mdata, 
    neighbors_key='wnn', 
    min_dist=0.3)

mu.write(snakemake.output.merged_multiome, compression='gzip')
