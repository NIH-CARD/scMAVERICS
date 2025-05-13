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

cell_type = snakemake.params.cell_type

adata = sc.read_h5ad(snakemake.input.celltype_atac)

adata.layers['counts'] = adata.X

# ATAC pseudobulk
pdata = dc.get_pseudobulk(
    adata,
    sample_col='sample_id',
    groups_col='cell_type',
    layer='counts',
    min_cells=1,
    min_counts=1
)

# 
pdata.layers['pcounts'] = pdata.X
sc.pp.normalize_total(pdata, target_sum=1e6, max_fraction = 0.000001, key_added='cpm', layer='pcounts', copy=False)

# 
covariate_df = pd.read_csv(snakemake.input.covariates)
covariate_df.index = covariate_df['Unnamed: 0']
del covariate_df['Unnamed: 0']
covariate_df = covariate_df.T
covariate_df['SampleID'] = covariate_df.index
covariates = ' + '.join(covariate_df.columns.to_list()[:-1])

# CSV pseudobulk for linear regression (includes covariates)
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[['sample_id']]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)
adata_df = pd.merge(left=covariate_df, right=adata_df, left_on='SampleID', right_on='sample_id')

#
adata_df.to_csv(snakemake.output.cell_specific_pseudo, index=False)

# 
peak_data = []
# 
for peak in pdata.var_names.to_list():
    #print(gene)
    # Create independent variable Series object, save as X
    _, X = dmatrices(f'Age ~ '+covariates, data=adata_df, return_type='dataframe')
    # Create dependent variable Series object, save as y
    y = adata_df[peak].values
    # Fit the model
    model = sm.OLS(y, X).fit(disp=0)
    # Add values to the dataframe
    peak_results = {
        'peak': peak,
        'slope': model.params['Age'],
        'p_value': model.pvalues['Age'],
        'standard_error': model.bse['Age'],
        'r_squared': model.rsquared,
        'adjusted_r_squared': model.rsquared_adj
    }
    peak_data.append(peak_results)

#
regression_df = pd.DataFrame(peak_data)

#
regression_df['p_value_bh'] = scipy.stats.false_discovery_control(regression_df['p_value'].fillna(1))
regression_df['-log10(p-value_bh)'] = -np.log10(regression_df['p_value_bh'])

# 
regression_df.to_csv(snakemake.output.cell_specific_regression, index=False)