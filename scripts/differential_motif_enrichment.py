import scanpy as sc
import numpy as np
import muon as mu
import pandas as pd
import scipy
import statsmodels.api as sm
import statsmodels.formula.api as smf

pchromvar sc.read_h5ad(snakemake.input.pseudobulk_chromvar)

chromvar.obs['Sample-celltype'] = chromvar.obs['Sample_ID'].astype(str) + '_' + chromvar.obs['celltype'].astype(str)

# Create Chromvar DataFrame
chromvar_df = chromvar.to_df()
chromvar_df['celltype'] = [x.split('_')[-1] for x in chromvar_df.index]
chromvar_df['Sample_ID'] = ['_'.join(x.split('_')[:-1]) for x in chromvar_df.index]

# Map diagnosis back onto DataFrame from RNA data
sample_diagnosis_dict = dict(zip(rna_adata_df['Sample_ID'].astype(str), rna_adata_df['Primary Diagnosis'].astype(str)))
chromvar_df['diagnosis'] = [sample_diagnosis_dict[x] for x in chromvar_df['Sample_ID']]

# Save list of TF motifs
TF_motif_names = chromvar_df.columns.to_list()
TF_motif_names = ['_'.join(x.replace("::", "..").replace("-", ".").split('.')) for x in TF_motif_names]
chromvar_df = chromvar_df.rename(columns=dict(zip(chromvar_df.columns.to_list(), TF_motif_names)))

chrom_rna_df = pd.merge(
    left = prna.to_df(),
    right = chromvar_df,
    left_index = True,
    right_index = True
)

disease_control = ['control', 'PD', 'LBD']
for i, condition_1 in enumerate(disease_control):
    for j, condition_2 in enumerate(disease_control):
        if i > j:
            comparison_combinations.append([condition_1, condition_2])

# List to append gene-motif links in 
gene_motif = []#= pd.DataFrame()
for celltype in cell_types:
    for comparisons in [['control', 'PD'], ['control', 'LBD'], ['PD', 'LBD']]:

        if comparisons[0] == 'control':
            comparison = comparisons[1]
        else:
            comparison = f'{comparisons[0]} vs. {comparison[1]}'
        
        # How many genes and DEMs to correlate
        condition_genes = sign_DEG_df[(sign_DEG_df['celltype'] == celltype) & (sign_DEG_df['diagnosis'] == comparison)]
        condition_motifs = sign_DEM_df[(sign_DEM_df['celltype'] == celltype) & (sign_DEM_df['comparison'] == comparison)]
        if condition_genes.shape[0] != 0 and condition_motifs.shape[0] != 0:
            print(celltype, comparison, condition_genes.shape[0], condition_motifs.shape[0])

            # For each gene-motif in both conditions, measure the Spearman correlation
            for diagnosis in comparisons:
                # Filter for cell type, diagnosis
                celltype_comparison_df = chrom_rna_df[(chrom_rna_df['celltype'] == celltype) & (chrom_rna_df['diagnosis'] == diagnosis)]
                for gene in condition_genes.id:
                    for motif in condition_motifs['TF motif']:
                        gene_express = celltype_comparison_df[gene]
                        motif_express = celltype_comparison_df[motif]
                        stat, p = scipy.stats.spearmanr(gene_express, motif_express)

                        gene_motif.append([celltype, diagnosis, gene, motif, stat, p])
gene_motif_df = pd.DataFrame(gene_motif, columns = ['celltype', 'diagnosis', 'gene', 'TF motif', 'Spearman rho', 'Spearman p-value'])