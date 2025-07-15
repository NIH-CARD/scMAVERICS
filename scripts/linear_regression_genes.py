import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import decoupler as dc
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from patsy import dmatrices

adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

sample_key = snakemake.params.sample_key
disease_param = snakemake.params.disease_param

# RNA pseudobulk
pdata = dc.get_pseudobulk(
    adata,
    sample_col=sample_key,
    groups_col='cell_type',
    layer='counts',
    mode='sum',
    min_cells=1,
    min_counts=1
)

pdata.layers['pcounts'] = pdata.X
sc.pp.normalize_total(pdata, target_sum=1e6, max_fraction = 0.001, key_added='cpm', layer='pcounts', copy=False)

# CSV pseudobulk for QTLs
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[[sample_key, 'cell_type']]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)

# CSV pseudobulk for linear regression (includes covariates)
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[snakemake.params.design_factors]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)

# Write out pseudobulk data
adata_df.to_csv(snakemake.output.rna_pseudobulk)

gene_data = []
for cell_type in adata_df['cell_type'].drop_duplicates():
    cohort_df = adata_df[(adata_df['cell_type'] == cell_type)]
    for gene in pdata.var_names.to_list():
        # Create independent variable Series object, save as X
        _, X = dmatrices(f'{disease_param} ~ Sex_numeric + PMI + Brain_bank_numeric', data=cohort_df, return_type='dataframe')
        # Create dependent variable Series object, save as y
        y = cohort_df[gene].values
        # Fit the model
        model = sm.OLS(y, X).fit(disp=0)
        # Add values to the dataframe
        gene_results = {
            'cell_type': cell_type,
            'gene': gene,
            'slope': model.params[disease_param],
            'p_value': model.pvalues[disease_param],
            'standard_error': model.bse[disease_param],
            'r_squared': model.rsquared,
            'adjusted_r_squared': model.rsquared_adj
        }
        gene_data.append(gene_results)

regression_df = pd.DataFrame(gene_data)
# Adjust p-value
for cell_type in adata_df['cell_type'].drop_duplicates():
    regression_df.loc[regression_df['cell_type'] == cell_type, 'p_value_bh'] = scipy.stats.false_discovery_control(regression_df.loc[regression_df['cell_type'] == cell_type, 'p_value'].fillna(1))
regression_df['-log10(p-value_bh)'] = -np.log10(regression_df['p_value_bh'])
regression_df.to_csv(snakemake.output.cell_gene_regression)