import muon as mu
import scanpy as sc


mdata = mu.MuData(
    {
        'rna': sc.read_h5ad(snakemake.input.merged_rna_anndata), # type: ignore
        'atac': sc.read_h5ad(snakemake.input.merged_atac_anndata) # type: ignore
    }
)



mu.pp.neighbors(mdata, key_added='wnn', metric='cosine', low_memory=True)

mu.tl.umap(mdata, neighbors_key='wnn', min_dist=0.3)


mdata.obsm['X_wnn_umap'] = mdata.obsm['X_umap']
mdata.mod['rna'].obsm['X_wnn_umap'] = mdata.obsm['X_wnn_umap']
mdata.mod['atac'].obsm['X_wnn_umap'] = mdata.obsm['X_wnn_umap']


mdata.write_h5mu(filename=snakemake.output.merged_multiome, compression='gzip') # type: ignore