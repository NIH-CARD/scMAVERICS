import muon as mu
import scanpy as sc

mdata = mu.MuData({"rna": sc.read_h5ad(snakemake.input.merged_rna_anndata), "atac": sc.read_h5ad(snakemake.input.merged_atac_anndata)})
mu.pp.intersect_obs(mdata)

mdata.write_h5mu(snakemake.output.merged_multiome)

