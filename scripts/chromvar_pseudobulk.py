import muon as mu
import scanpy as sc

# Read in MuData
mdata = mu.read_h5mu(snakemake.input.merged_multiome)

# Separate out the chromvar modality
chromvar  = mdata.mod['chromvar']

# Add possibly missing data
for parameter in [snakemake.params.sample_key, 'celltype']:
    paramdict = mdata.obs[parameter].to_dict()
    chromvar.obs[parameter] = [paramdict[x] for x in chromvar.obs_names]

chromvar.obs['Sample-celltype'] = chromvar.obs[snakemake.params.sample_key].astype(str) + '_' + chromvar.obs[snakemake.params.separating_cluster].astype(str)
pchromvar = sc.get.aggregate(
    chromvar, 
    by='Sample-celltype', 
    func = 'mean'
    )
pchromvar.X = pchromvar.layers['mean']

pchromvar.write_h5ad(snakemake.output.pseudobulk_chromvar, compression='gzip')
