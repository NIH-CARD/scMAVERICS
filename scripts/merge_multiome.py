import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import scipy.sparse as sp

# Read in RNA dataset
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Read in ATAC dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Save the anndata object
mdata = mu.MuData({"rna": rna, "atac": atac})

mdata.write(snakemake.output.multiome_object, compression ='gzip')
