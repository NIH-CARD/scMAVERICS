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

# RNA pseudobulk
pdata = dc.get_pseudobulk(
    adata,
    sample_col='SampleID',
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
sample_cell = pdata.obs[['SampleID', 'cell_type']]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)

for cell_type in ['Astro', 'ExN', 'InN', 'MG', 'OPC', 'Oligo', 'VC']:
    cell_adata_df = adata_df[adata_df['cell_type'] == cell_type]
    cell_adata_df.index = cell_adata_df['SampleID']
    del cell_adata_df['SampleID']
    del cell_adata_df['cell_type']
    cell_adata_df.to_csv(f'/data/CARD_singlecell/PFC_atlas/data/celltypes/{cell_type}/pseudobulk_rna.csv')

# CSV pseudobulk for linear regression (includes covariates)
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[['SampleID', 'cell_type', 'PMI', 'Sex', 'Brain_bank', 'Age']]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)

adata_df['PMI'] = adata_df['PMI'].astype(float)
adata_df['Age'] = adata_df['Age'].astype(float)

# Write out pseudobulk data
adata_df.to_csv(snakemake.output.rna_pseudobulk)

sex2num = dict(zip(adata_df['Sex'].drop_duplicates().to_list(), range(len(adata_df['Sex'].drop_duplicates().to_list()))))
adata_df['Sex_numeric'] = [sex2num[x] for x in adata_df['Sex']]

bbank2num = dict(zip(adata_df['Brain_bank'].drop_duplicates().to_list(), range(len(adata_df['Brain_bank'].drop_duplicates().to_list()))))
adata_df['Brain_bank_numeric'] = [bbank2num[x] for x in adata_df['Brain_bank']]

gene_data = []
for cell_type in adata_df['cell_type'].drop_duplicates():
    cohort_df = adata_df[(adata_df['cell_type'] == cell_type)]
    for gene in pdata.var_names.to_list():
        #print(gene)
        # Create independent variable Series object, save as X
        _, X = dmatrices(f'Age ~ Age + Sex_numeric + PMI + Brain_bank_numeric', data=cohort_df, return_type='dataframe')
        # Create dependent variable Series object, save as y
        y = cohort_df[gene].values
        # Fit the model
        model = sm.OLS(y, X).fit(disp=0)
        # Add values to the dataframe
        gene_results = {
            'cell_type': cell_type,
            'gene': gene,
            'slope': model.params['Age'],
            'p_value': model.pvalues['Age'],
            'standard_error': model.bse['Age'],
            'r_squared': model.rsquared,
            'adjusted_r_squared': model.rsquared_adj
        }
        gene_data.append(gene_results)

regression_df = pd.DataFrame(gene_data)
# Adjust p-value
for cell_type in adata_df['cell_type'].drop_duplicates():
    regression_df.loc[regression_df['cell_type'] == cell_type, 'p_value_bh'] = scipy.stats.false_discovery_control(regression_df.loc[regression_df['cell_type'] == cell_type, 'p_value'].fillna(1))
regression_df['-log10(p-vlue_bh)'] = -np.log10(regression_df['p_value_bh'])
regression_df.to_csv(snakemake.output.cell_specific_regression)