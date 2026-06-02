import muon as mu
import scanpy as sc

chromvar  = mu.read(snakemake.input.merged_multiome + "/mod/chromvar")

chromvar.obs['Sample-celltype'] = chromvar.obs[snakemake.params.sample_key].astype(str) + '_' + chromvar.obs[snakemake.params.separating_cluster].astype(str)
pchromvar = sc.get.aggregate(
    chromvar, 
    by='Sample-celltype', 
    func = 'mean'
    )
pchromvar.X = pchromvar.layers['mean']

pchromvar.write_h5ad(snakemake.output.pseudobulk_chromvar, compression='gzip')