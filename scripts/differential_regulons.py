DEG_df = pd.read_csv('../../data/DGE_Dreampy_results.csv')
sign_DEG_df = DEG_df[
    (DEG_df['coef'].isin(['Primary Diagnosis_control', 'Primary Diagnosis_PD'])) &
    (DEG_df['adj_p_val'] < 0.05) &
    (abs(DEG_df['logfc']) > 0.5)]
sign_DEG_df


for celltype in cell_types:
    for diagnosis in ['control', 'PD' ,'LBD']:
        try:
            celltype_gene_motif = pd.read_csv(work_dir + f'/data/celltypes/{celltype}/{celltype}_{diagnosis}_motif_gene_connect.csv')

celltype_gene_peak_linkage = pd.read_csv(work_dir + f'/data/celltypes/{celltype}/{celltype}_promoter_coaccessibility.csv')

celltype_gene_peak_linkage_df = celltype_gene_peak_linkage_df[['gene name', 'diagnosis', 'score', 'overlap promoter peak', 'overlap enhancer peak', 'gene-peak link stat', 'gene-peak link p-value', 'celltype']].drop_duplicates()