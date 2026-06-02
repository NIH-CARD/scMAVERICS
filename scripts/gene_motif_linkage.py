import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr,
from statsmodels.stats.multitest import multipletests

# Import the celltype of interest
celltype_param = snakemake.params.celltype_param
celltype = snakemake.params.celltype
# Import the diagnosis of interest, and the parameter that the DataFrame is sliced from
diagnosis_param = snakemake.params.diagnosis_param
diagnosis = snakemake.params.diagnosis

# Import pseudobulked chromvar TF motif data
pchromvar = sc.read_h5ad(snakemake.input.pseudobulk_chromvar)
# Make sure the imported raw layer is the mean aggregated signal
pchromvar.X = pchromvar.layers['mean']

# Import pseudobulked rna data
prna = sc.read_h5ad(snakemake.input.pseudo_rna)
# Make sure the imported raw data is also the correct, log-normalized signal
prna.X = prna.layers['X']

# Convert the Anndata object into a DataFrame
pchromvar_df = pchromvar.to_df()

# Save the motif names to a list, convert them be separated by _
TF_motif_names = pchromvar_df.columns.to_list()
TF_motif_names = ['_'.join(x.replace("::", "..").replace("-", ".").split('.')) for x in TF_motif_names]
# Rename the columns with the corrected TF motif signals
pchromvar_df = pchromvar_df.rename(columns=dict(zip(pchromvar_df.columns.to_list(), TF_motif_names)))

# Combine the RNA and TF motif signals
prna_df = prna.to_df()
gene_names = prna_df.columns.to_list()
gene_motif_df = pd.merge(
    left = prna_df,
    right = pchromvar_df,
    left_index = True,
    right_index = True
)
gene_motif_df = pd.merge(
    left = gene_motif_df,
    right = prna.obs[['celltype', diagnosis_param]],
    left_index = True,
    right_index = True
)

diagnosis_gene_motif_df = pd.DataFrame()
celltype_matrix_df = gene_motif_df[(gene_motif_df['celltype'] == celltype) & (gene_motif_df['Primary Diagnosis'] == diagnosis)]
for motif in TF_motif_names:
    
    res = celltype_matrix_df[gene_names].apply(lambda gene_vec: pearsonr(gene_vec, celltype_matrix_df[motif]))
    correlation_df = res.iloc[[0]].melt(var_name ='gene', value_name= 'Pearson r')
    p_value_df = res.iloc[[1]].melt(var_name ='gene', value_name= 'Pearson p-value')

    temp_df = pd.merge(
        left = correlation_df,
        right = p_value_df,
        left_on = 'gene',
        right_on = 'gene'
    )
    temp_df['motif'] = motif
    temp_df['celltype'] = celltype
    temp_df['diagnosis'] = diagnosis
    diagnosis_gene_motif_df = pd.concat([diagnosis_gene_motif_df, temp_df])

# Correct for false discovery with Bonferroni-Holm method, by diagnosis
diagnosis_gene_motif_df['adj. p-value'] = multipletests(
    pvals = diagnosis_gene_motif_df['Pearson p-value'],
    alpha = 0.05,
    method = 'holm')[1]
# Add diagnosis based group to sample
diagnosis_gene_motif_df['-log10(adj. p-value)'] = -np.log10(diagnosis_gene_motif_df['adj. p-value'])

# Export gene-motif linkages
diagnosis_gene_motif_df.to_csv(snakemake.output.gene_motif_links, index = False)
