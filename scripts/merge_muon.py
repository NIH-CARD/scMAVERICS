import muon as mu
import scanpy as sc

mdata = mu.MuData({"RNA": sc.read_h5ad(snakemake.input.merged_rna_anndata), "ATAC": sc.read_h5ad(snakemake.input.merged_atac_anndata)})

mdata.write_h5mu(snakemake.output.merged_multiome)
