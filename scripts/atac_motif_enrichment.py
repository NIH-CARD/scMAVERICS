
chromvar  = mu.read("../../atlas/multiome_chromvar_atlas.h5mu/mod/chromvar")
chromvar

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

DEM_df[(DEM_df['adj. p-value'] < 0.05) & (abs(DEM_df['log2FC']) > 0.5)].groupby(['celltype', 'comparison']).count()
