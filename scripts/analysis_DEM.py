import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import pysam
import collections
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import multitest

# parameters to add
sample_key = snakemake.params.sample_key
min_nuclei = snakemake.params.min_nuclei


# Define parameters
separating_cluster = snakemake.params.separating_cluster
diagnosis_control = diagnoses
control_param = snakemake.params.control
comparison_param = snakemake.params.diagnosis_param
random_cov_var = snakemake.params.random_cov_var
fixed_cont_var = snakemake.params.fixed_cont_var
fixed_cate_var = snakemake.params.fixed_cate_var

# Read in pseudobulked chromvar data
pchromvar = sc.read_h5ad(snakemake.input.pseudobulked_chromvar)
chromvar_df = pchromvar.to_df()

TF_motif_names = chromvar_df.columns.to_list()
TF_motif_names = ['_'.join(x.replace("::", "..").replace("-", ".").split('.')) for x in TF_motif_names]
chromvar_df = chromvar_df.rename(columns=dict(zip(chromvar_df.columns.to_list(), TF_motif_names)))
chromvar_df[sample_key] = ['_'.join(x.split('_')[:-1]) for x in chromvar_df.index]
chromvar_df[separating_cluster] = [x.split('_')[-1] for x in chromvar_df.index]

# Read in metadata and merge
metadata_df = pd.read_csv(snakemake.input.metadata)

chromvar_df = pd.merge(
    left = chromvar_df,
    right = metadata_df,
    left_on = sample_key,
    right_on = sample_key
)

DEM_df = pd.DataFrame()
for celltype in chromvar_df.celltype.drop_duplicates():
    print(celltype)
    for comparison in comparison_combinations:

        if comparison[0] == control_param:
            disease_name = comparison[1]
        else:
            disease_name = f'{comparison[1]} vs. {comparison[0]}'
        celltype_comparison_chromvar_df = chromvar_df[(chromvar_df[separating_cluster] == celltype) & (chromvar_df[comparison_param].isin(comparison))]
        
        # Don't run if lack values
        if celltype_comparison_chromvar_df.groupby(comparison_param).count().shape[0] > 1:

            # First add continuous variables
            FORMULA = f' ~ {comparison_param}'
            for cont_var in fixed_cont_var:
                FORMULA += f' + {cont_var}'

            # First add continuous variables
            for cont_var in fixed_cont_var:
                FORMULA += f' + {cont_var}'

            # Random covariates conversion sub-script
            cov_column = []
            if len(random_cov_var) != 0:
                random_list = [celltype_comparison_chromvar_df[column].astype(str).to_list() for column in random_cov_var]
                for i in range(len(random_list[0])):
                    new_line_cov = []
                    for j in range(len(random_cov_var)):
                        new_line_cov.append(random_list[j][i])
                    cov_column.append('_'.join(new_line_cov))
                celltype_comparison_chromvar_df['random_cov'] = cov_column
                celltype_comparison_chromvar_df['random_cov'] = celltype_comparison_chromvar_df['random_cov'].astype('category')
            
            celltype_motif_slope_list = []
            for TF_motif in TF_motif_names:
                TF_FORMULA = TF_motif + FORMULA
                if len(random_cov_var) != 0:
                    ccc_model = smf.mixedlm(TF_FORMULA, celltype_comparison_chromvar_df, groups = celltype_comparison_chromvar_df['random_cov'])
                    mdf = ccc_model.fit(method=["lbfgs"])
                else:
                    ccc_model = smf.glm(TF_FORMULA, celltype_comparison_chromvar_df)
                    mdf = ccc_model.fit()
                result_param = mdf.params[mdf.params.index.str.contains(comparison_param)].index[0]
                celltype_motif_slope_list.append([f"{TF_motif}", mdf.params[result_param], mdf.pvalues[result_param], mdf.params.Intercept])

            celltype_motif_slope_df = pd.DataFrame(celltype_motif_slope_list, columns = ['TF motif', 'log2FC', 'p-value', 'intercept'])
            celltype_motif_slope_df['adj. p-value'] = multitest.multipletests(pvals = celltype_motif_slope_df['p-value'], alpha=0.01, method = 'holm')[1]
            celltype_motif_slope_df['-log10(adj. p-value)'] = -np.log10(celltype_motif_slope_df['adj. p-value'])

            # Add celltype and condition specific parameters
            celltype_motif_slope_df[separating_cluster] = celltype
            celltype_motif_slope_df['comparison'] = disease_name

            DEM_df = pd.concat([DEM_df, celltype_motif_slope_df])

DEM_df.to_csv(snakemake.output.diff_enrich_motif, compression = 'gzip')