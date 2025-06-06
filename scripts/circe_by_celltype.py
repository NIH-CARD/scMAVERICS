import numpy as np
import circe as ci
import scanpy as sc
import scipy as sp

adata = sc.read_h5ad(snakemake.input.celltype_atac)
adata.var['start'] = adata.var['start'].astype(int)
adata.var['end'] = adata.var['end'].astype(int)
adata.var['chr'] = adata.var['chromosome']

adata = ci.metacells.compute_metacells(
    adata, 
    method='sum', 
    k=25, 
)

# Get co-accessibility scores
final_score = ci.sliding_graphical_lasso(
    adata,
    n_samples=50,
    n_samples_maxtry=100,
    max_alpha_iteration=500,
    verbose=True
)

adata.varp['atac_network'] = final_score

adata.write_h5ad(snakemake.output.celltype_atac, compression='gzip')

# Compute the co-accessibility network
atac = ci.add_region_infos(adata, sep=(':', '-'))
ci.compute_atac_network(atac)

# Extract the network and find CCANs modules
circe_network = ci.extract_atac_links(atac)

circe_network['Peak1'] = circe_network['Peak1'].str.replace(':', '_').str.replace('-', '_')
circe_network['Peak2'] = circe_network['Peak2'].str.replace(':', '_').str.replace('-', '_')

circe_network.to_csv(snakemake.output.circe_network, index=None)