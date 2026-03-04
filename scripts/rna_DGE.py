import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Open rna
adata = sc.read_h5ad(snakemake.input.rna_anndata)

# Read in parameters
cell_type = snakemake.params.cell_type
disease_name = snakemake.params.disease
control_name = snakemake.params.control
disease_param = snakemake.params.disease_param

# Subset to cell type
adata = adata[adata.obs[snakemake.params.separating_cluster] == cell_type].copy()

# Get pseudo-bulk profile
pdata = dc.pp.pseudobulk(
    SN_adata,
    sample_col=snakemake.params.sample_key,
    groups_col=snakemake.params.separating_cluster,
    layer='counts',
    mode='sum'
)

dc.pp.filter_samples(pdata, min_cells=10, min_counts=1000)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

dc.pl.filter_by_expr(
    adata=pdata,
    group=disease_param,
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
)
dc.pl.filter_by_prop(
    adata=pdata,
    min_prop=0.1,
    min_smpls=2,
)

# Abbreviate diagnosis to avoid space syntax error
pdata.obs['comparison'] = pdata.obs[disease_param]

# Determine the number of cpus to use
inference = DefaultInference(n_cpus=64)

# Design the differential expression analysis with covariates
dds = DeseqDataSet(
    adata=pdata,
    design_factors=snakemake.params.design_factors + ['comparison'],
    refit_cooks=True,
    inference=inference,
)

# Error handling for singular matrixes (they don't vary across co-variates)
try:
    # Compute log-fold changes
    dds.deseq2()

    # Extract contrast between control and disease states
    stat_res = DeseqStats(
        dds,
        contrast=['comparison', disease_name, control_name],
        inference=inference,
    )

    # Compute Wald test
    stat_res.summary()

    # Extract results
    DGE_results_df = stat_res.results_df
    DGE_results_df[snakemake.params.separating_cluster] = cell_type
    DGE_results_df['-log10_padj'] = -np.log10(DGE_results_df['padj'])
    DGE_results_df.to_csv(snakemake.output.output_DGE_data)

    # Plot 
    dc.plot_volcano_df(
        DGE_results_df,
        x='log2FoldChange',
        y='padj',
        top=20,
        lFCs_thr=1,
        sign_thr=1e-2,
        figsize=(4, 4)
    )
    plt.title(f'{control_name} vs. {disease_name} in {cell_type}')
    plt.tight_layout()
    plt.savefig(snakemake.output.output_figure, dpi=300)
except np.linalg.LinAlgError as err:
    print('\n\nSingular matrix! No output\n\n')
    # Output nothing if 
    pd.DataFrame().to_csv(snakemake.output.output_DGE_data)
    plt.plot()
    plt.savefig(snakemake.output.output_figure, dpi=300)

    