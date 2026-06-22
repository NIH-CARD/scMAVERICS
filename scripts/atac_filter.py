import scanpy as sc
import snapatac2 as snap

# Load in data 

# Important that ATAC is read in using snapATAC2
atac = snap.read(snakemake.input.atac_anndata, backed=None) 
# Filter cells with less than 1000 counts
atac = atac[((atac.obs['n_fragment'] > snakemake.params.min_peak_counts) & (atac.obs['tsse'] > snakemake.params.min_tsse))].copy()

atac.write(snakemake.output.atac_anndata, compression='gzip')