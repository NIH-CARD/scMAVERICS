import scanpy as sc
import snapatac2 as snap

# Load in data 

# Important that ATAC is read in using snapATAC2
atac = snap.read(snakemake.input.atac_anndata, backed=None) 
# RNA can be read in the normal Scanpy way
rna = sc.read_h5ad(snakemake.input.rna_anndata)

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
atac.write(snakemake.output.atac_anndata, compression='gzip')
