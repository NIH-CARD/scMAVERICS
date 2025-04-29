from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Read in annotated ATAC data
atac = sc.read_h5ad('/data/CARD_singlecell/SN_atlas/atlas/04_modeled_anndata_atac.h5ad')

# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    atac,
    sample_col='sample_id',
    groups_col='cell_type',
    mode='sum',
    min_cells=10,
    min_counts=10
    )

# CSV pseudobulk
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[['sample_id']]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)
adata_df.to_csv('/data/CARD_AUX/SN_PFC_MULTIOME/SN/SN_RNA_atlas/snATAC_pseudobulk.csv', index=False)

# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()

# Abreviate diagnosis
pdata.obs['diagnosis'] = pdata.obs['Primary Diagnosis']

# Select Astrocytes profiles
astrocytes = pdata[pdata.obs['cell_type'] == 'Astro'].copy()
astrocyte_genes = dc.filter_by_expr(astrocytes, group='diagnosis', min_count=10, min_total_count=15)
astrocytes = astrocytes[:, astrocyte_genes].copy()

inference = DefaultInference(n_cpus=32)

dds = DeseqDataSet(
    adata=astrocytes,
    design_factors=['normalage', 'diagnosis'],
    inference=inference,
)

# Compute LFCs
dds.deseq2()

# Extract contrast between normal and DLB
stat_res = DeseqStats(
    dds,
    contrast=["diagnosis", 'DLB', 'control'],
    inference=inference,
)

# Compute Wald test
stat_res.summary()

# Extract results
astrocyte_DLB_results_df = stat_res.results_df

dc.plot_volcano_df(
    astrocyte_DLB_results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    lFCs_thr=1,
    sign_thr=1e-2,
    figsize=(4, 4)
)
plt.title('Control vs. DLB in Astrocytes')
plt.tight_layout()
plt.savefig('/data/CARD_singlecell/SN_atlas/figures/Astro/DLB_DAR_volcano.svg', dpi=300)

# Write out dataframe containing DARs
