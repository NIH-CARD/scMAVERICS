import muon as mu
import scanpy as sc

# Read in data
mdata = mu.read_h5mu(snakemake.input.merged_multiome)

# Split out only the celltype AnnData, and reduce the size of the object
mdata = mdata[mdata['RNA'].obs['cell_type'] == snakemake.params.cell_type].copy()

# Split out the RNA AnnData
rna = mdata['RNA'].copy()

# Split out the ATAC AnnData
atac = mdata['ATAC'].copy()

# Write the rna data out
rna.write_h5ad(snakemake.output.celltype_rna, compression='gzip')

# Write the atac data out
atac.write_h5ad(snakemake.output.celltype_atac, compression='gzip')
