import scanpy as sc

# Read in the annotated gene expression dataset
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Read in the overlapping chromatin accessibility dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Create dictionary of nuclei ID to cell type annotation from the
# gene expression dataset
cell_annot = adata.obs['celltype'].to_dict()

# Create list of cell type annotations that matches the list of atac indexes
atac.obs[snakemake.params.pseudobulk_param] =  [cell_annot[x] for x in atac.obs_names]

# Export atac datafile
atac.write_h5ad(
    filename=snakemake.output.merged_atac_anndata, 
    compression='gzip'
)