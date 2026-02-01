import pandas as pd
import numpy as np
import scanpy as sc
import liana as li

# Open the AnnData object
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)
# Subset to disease state
PD_adata = adata[adata.obs[snakemake.params.disease_param].isin([snakemake.params.disease])].copy()

# Open previous DGE calculation
subtype_DGE_df = pd.read_csv(snakemake.input.merged_DGE_data)

disease_subtype_DGE_df = subtype_DGE_df[subtype_DGE_df[snakemake.params.disease_param] == snakemake.params.disease]
disease_subtypes = disease_subtype_DGE_dfp[snakemake.params.sep_param].drop_duplicates().to_list()

# Dictionary to store gene DGE data, separated by cell separating factor
subtype_dea = {}
for subtype in disease_subtypes:
    subtype_df = disease_subtype_DGE_df[disease_subtype_DGE_df[snakemake.params.sep_param] == subtype]
    subtype_df = subtype_df.set_index('gene')
    subtype_dea[subtype] = subtype_df
dea_df = pd.concat(subtype_dea)
del dea_df[snakemake.params.sep_param]
dea_df = dea_df.reset_index().rename(columns={'level_0': snakemake.params.sep_param}).set_index('gene')

lr_res = li.multi.df_to_lr(
    PD_adata,
    dea_df=dea_df,
    resource_name='consensus', # NOTE: uses HUMAN gene symbols!
    expr_prop=0.1, # calculated for adata as passed - used to filter interactions
    groupby=snakemake.params.sep_param,
    stat_keys=['stat', 'pvalue', 'padj'],
    use_raw=False,
    complex_col='stat', # NOTE: we use the Wald Stat to deal with complexes
    verbose=True,
    return_all_lrs=False,
    )

lr_res.to_csv(snakemake.output.differential_cell_cell_communication_data, index=False)