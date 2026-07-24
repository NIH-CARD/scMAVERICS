import scanpy as sc
import numpy as np
import muon as mu
import pandas as pd
import scipy
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Load RNA data
rna_adata_df = pd.read_csv('../../data/rna_cell_annot.csv')
rna_adata_df.index = rna_adata_df['Unnamed: 0']
del rna_adata_df['Unnamed: 0']
rna_adata_df

chromvar  = mu.read("../../atlas/multiome_chromvar_atlas.h5mu/mod/chromvar")


chromvar.obs = pd.merge(
    left = chromvar.obs,
    right = rna_adata_df,
    left_index = True,
    right_index = True
)

# Merged sample-celltype name inputs
sample_key = snakemake.parms.sample_key
celltype = snakemake.params.separating_cluster

chromvar.obs['Sample-celltype'] = chromvar.obs[sample_key].astype(str) + '_' + chromvar.obs[celltype].astype(str)

pchromvar = sc.get.aggregate(
    chromvar, 
    by='Sample-celltype', 
    func = 'mean'
    )
pchromvar.X = pchromvar.layers['mean']

# Create Chromvar DataFrame
chromvar_df = pchromvar.to_df()
chromvar_df['celltype'] = [x.split('_')[-1] for x in chromvar_df.index]
chromvar_df['Sample_ID'] = ['_'.join(x.split('_')[:-1]) for x in chromvar_df.index]

# Map diagnosis back onto DataFrame from RNA data
sample_diagnosis_dict = dict(zip(rna_adata_df['Sample_ID'].astype(str), rna_adata_df['Primary Diagnosis'].astype(str)))
chromvar_df['diagnosis'] = [sample_diagnosis_dict[x] for x in chromvar_df['Sample_ID']]

# Save list of TF motifs
TF_motif_names = chromvar_df.columns.to_list()
TF_motif_names = ['_'.join(x.replace("::", "..").replace("-", ".").split('.')) for x in TF_motif_names]
chromvar_df = chromvar_df.rename(columns=dict(zip(chromvar_df.columns.to_list(), TF_motif_names)))

# Merge data 
chrom_rna_df = pd.merge(
    left = prna.to_df(),
    right = chromvar_df,
    left_index = True,
    right_index = True
)

DEM_df = pd.DataFrame()
for celltype in chromvar_df.celltype.drop_duplicates():
    for comparison in [['control', 'PD'], ['control', 'LBD'], ['PD', 'LBD']]:

        if comparison[0] == 'control':
            disease_name = comparison[1]
        else:
            disease_name = f'{comparison[1]} vs. {comparison[0]}'
        celltype_comparison_chromvar_df = chromvar_df[(chromvar_df['celltype'] == celltype) & (chromvar_df['diagnosis'].isin(comparison))]
        celltype_comparison_chromvar_df['batch_bank'] = celltype_comparison_chromvar_df['Use_batch'].astype(str) + '-' + celltype_comparison_chromvar_df['Brain_bank'].astype(str)
        celltype_comparison_chromvar_df['batch_bank'] = celltype_comparison_chromvar_df['batch_bank'].astype('category')

        # One-hot encode
        celltype_comparison_chromvar_df['diagnosis_onehot'] = [1 if x == comparison[0] else 0 for x in celltype_comparison_chromvar_df.diagnosis]
        celltype_comparison_chromvar_df['Sex_onehot'] = [1 if x == 'Male' else 0 for x in celltype_comparison_chromvar_df.Sex]
        celltype_motif_slope_list = []
        for TF_motif in TF_motif_names:
            ccc_model = smf.mixedlm(f"{TF_motif} ~ diagnosis_onehot + Age + Sex_onehot", celltype_comparison_chromvar_df, groups = celltype_comparison_chromvar_df['batch_bank'])
            mdf = ccc_model.fit(method=["lbfgs"])
            celltype_motif_slope_list.append([f"{TF_motif}", mdf.params.diagnosis_onehot, mdf.pvalues.diagnosis_onehot, mdf.params.Intercept])


            
        celltype_motif_slope_list

        celltype_motif_slope_df = pd.DataFrame(celltype_motif_slope_list, columns = ['TF motif', 'log2FC', 'p-value', 'intercept'])
        celltype_motif_slope_df['adj. p-value'] = multitest.multipletests(pvals = celltype_motif_slope_df['p-value'], alpha=0.01, method = 'holm')[1]
        celltype_motif_slope_df['-log10(adj. p-value)'] = -np.log10(celltype_motif_slope_df['adj. p-value'])

        # Add celltype and condition specific parameters
        celltype_motif_slope_df['celltype'] = celltype
        celltype_motif_slope_df['comparison'] = disease_name

        DEM_df = pd.concat([DEM_df, celltype_motif_slope_df])
DEM_df.to_csv(snakemake.output.diff_enrich_motif)