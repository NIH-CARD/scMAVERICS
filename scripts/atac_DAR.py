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
pdata = dc.pp.pseudobulk(
    celltype_atac,
    sample_col='sample_id',
    groups_col=snakemake.params.separating_cluster,
    mode='sum'
)

# Filter out samples with low number of cells or counts by each param and covariate
dc.pp.filter_samples(pdata, min_cells=10, min_counts=100)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

dc.pp.filter_by_expr(
    adata=pdata,
    group="disease",
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
)
dc.pp.filter_by_prop(
    adata=pdata,
    min_prop=0.1,
    min_smpls=2,
)

# Include inference
inference = DefaultInference(n_cpus=snakemake.threads)

# Design the differential expression analysis with covariates
dds = DeseqDataSet(
        adata=subtype_pdata,
        design_factors=[disease_param] + snakemake.params.design_covariates,
        refit_cooks=True,
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
deseq2_results_df['peak'] = deseq2_results_df.index.to_list()
deseq2_results_df.to_csv(snakemake.output.output_DAR_data, index=False)

# Plot
dc.pl.plot_volcano_df(
    deseq2_results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    thr_stat=1,
    thr_sign=1e-2,
    figsize=(4, 4)
)
plt.title(f'{control_name} vs. {disease_name} in {cell_type}')
plt.tight_layout()
plt.savefig(snakemake.output.output_figure, dpi=300)
