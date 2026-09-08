# Code for cell-cell communication
import liana as li
import pandas as pd
import numpy as np
import scanpy as sc

adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

control_adata = adata[adata.obs[snakemake.params.disease_param] == snakemake.params.control]

# Run rank_aggregate
li.mt.cellphonedb(
    control_adata,
    groupby='celltype',
    #resource=cellphone_resource,
    resource_name='consensus',
    expr_prop=0.1,
    verbose=True,
    key_added='li_res',
    n_jobs=snakemake.threads
    )

li.mt.rank_aggregate(control_adata, 
                     groupby='celltype',
                     resource_name='consensus',
                     expr_prop=0.1,
                     verbose=True)

cellphone_df = control_adata.uns['li_res'].copy()
cellphone_df.to_csv(snakemake.output.cell_cell_communication_data, index=None)