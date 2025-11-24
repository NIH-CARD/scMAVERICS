import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Read in annotated ATAC data
atac = sc.read_h5ad(snakemake.input.atac_anndata)

# Read in parameters
cell_type = snakemake.params.cell_type
disease_name = snakemake.params.disease
control_name = snakemake.params.control
disease_param = snakemake.params.disease_param

# Subset to cell type
atac = atac[atac.obs[snakemake.params.separating_cluster] == cell_type].copy()

# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    atac,
    sample_col='sample_id',
    groups_col=disease_param,
    mode='sum',
    min_cells=10,
    min_counts=10
    )

# CSV pseudobulk
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[['sample_id']]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df.to_csv(snakemake.output.cell_specific_pseudo, index=False)

# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()

# Abbreviate diagnosis to avoid space syntax error
pdata.obs['comparison'] = pdata.obs[disease_param]

# Select gene specific profiles
pdata_genes = dc.filter_by_expr(
    pdata, 
    group='comparison', 
    min_count=10, 
    min_total_count=15
    )

# Subset valuable genes
pdata = pdata[:, pdata_genes].copy()

# Include inference
inference = DefaultInference(n_cpus=64)

# Design the differential expression analysis with covariates
dds = DeseqDataSet(
    adata=pdata,
    design_factors=snakemake.params.design_factors + ['comparison'],
    inference=inference,
)

# Compute LFCs
dds.deseq2()

# Extract contrast between control and disease state
stat_res = DeseqStats(
    dds,
    contrast=['comparison', disease_name, control_name],
    inference=inference,
)

# Compute Wald test
stat_res.summary()

# Extract results
deseq2_results_df = stat_res.results_df
deseq2_results_df['-log10_padj'] = -np.log10(deseq2_results_df['padj'])
deseq2_results_df.to_csv(snakemake.output.output_DAR_data)

# Plot
dc.plot_volcano_df(
    deseq2_results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    lFCs_thr=1,
    sign_thr=1e-2,
    figsize=(4, 4)
)
plt.title(f'Control vs. {disease_name} in {cell_type}')
plt.tight_layout()
plt.savefig(snakemake.output.output_figure, dpi=300)
