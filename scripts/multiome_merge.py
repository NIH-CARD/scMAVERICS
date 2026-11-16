import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import scipy.sparse as sp
ad.settings.allow_write_nullable_strings = True


# Read in RNA dataset
rna = sc.read_h5ad(snakemake.input.merged_rna_anndata)
rna = rna[(rna.obs['celltype'] != 'NA') & ~(rna.obs['CARD_ID'] == 'PCA-008')].copy()

# Read in ATAC dataset
atac = sc.read_h5ad(snakemake.input.merged_atac_anndata)
atac = atac[~(atac.obs['sample'] == 'PCA-008') & (atac.obs['n_fragment'] > 4000)].copy()

mdata = mu.MuData({"rna": rna, "atac": atac})

mdata.obs[snakemake.params.sample_key] = ['-1'.join(x.split('-1_')[1:]) for x in mdata.obs_names]

mdata.write(snakemake.output.multiome_object)
