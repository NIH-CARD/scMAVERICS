import scanpy as sc

# Load in data 

# Important that ATAC is read in using snapATAC2
atac = sc.read_h5ad(snakemake.input.atac_anndata) 
# RNA can be read in the normal Scanpy way
rna = sc.read_h5ad(snakemake.input.rna_anndata)
# Filter cells with less than the input number of counts
atac = atac[atac.obs['n_fragment'] > snakemake.params.min_peak_counts].copy()
# Filter out based on minimum number of TSSE
atac = atac[atac.obs['tsse'] > snakemake.params.min_tsse].copy()

# Add the consolidated cell-barcode 'atlas_identifier'
rna_barcodes = rna.obs_names
atac_barcodes = atac.obs_names

# Subset RNA and ATAC objects based on the overlap of values
rna = rna[rna.obs_names.isin(atac_barcodes)].copy()
atac = atac[atac.obs_names.isin(rna_barcodes)].copy()

rna.write_h5ad(filename=snakemake.output.rna_anndata, compression='gzip')
atac.write_h5ad(filename=snakemake.output.atac_anndata, compression='gzip')
