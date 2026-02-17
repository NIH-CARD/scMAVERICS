import scanpy as sc
import snapatac2 as snap

# Load in data 

# Important that ATAC is read in using snapATAC2
atac = sc.read_h5ad(snakemake.input.atac_anndata) 
# RNA can be read in the normal Scanpy way
rna = sc.read_h5ad(snakemake.input.rna_anndata)
# Filter cells with less than the input number of counts
atac = atac[atac.obs['n_fragment'] > snakemake.params.min_peak_counts].copy()
# Filter out based on minimum number of TSSE
atac = atac[atac.obs['tsse'] > snakemake.params.min_tsse].copy()

# List of RNA barcodes
rna_in_atac = [x for x in rna.obs_names if x in atac.obs_names]

# Subset the ATAC file with only the barcoes in RNA list
atac_subset = atac.subset(
    obs_indices=rna_in_atac,
    out=snakemake.output.atac_anndata,
    inplace=False
    )

# Filter the cells with the same barcodes as the RNA sample
rna = rna[rna.obs.index.isin(atac.obs_names)].copy()
rna.write_h5ad(filename=snakemake.output.rna_anndata, compression='gzip')

# Close both unfiltered and filtered datasets
atac_subset.close()
atac.close()