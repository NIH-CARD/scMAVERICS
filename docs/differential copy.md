# scMAVERICS - differential analysis




```python
# Load packages

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import decoupler as dc
import dreampy as dp
import numpy as np
import pandas as pd
import gseapy as gp
import greatpy as great
```

## Pseudobulking




```python

```

### RNA pseudobulking




```python
# Open rna
adata = sc.read_h5ad(snakemake.input.rna_anndata)
```


```python
# Get pseudo-bulk profile, summing all raw counts by cell type and sample
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=snakemake.params.sample_key,
    groups_col=snakemake.params.separating_cluster,
    layer='counts',
    mode='sum'
)
```


```python


dc.pp.filter_samples(pdata, min_cells=snakemake.params.min_cells)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()
```


```python
# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata)
sc.tl.pca(pdata)
```


```python
# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

dc.pp.filter_by_expr(
    adata=pdata,
    group="celltype",
    min_count=10,
    min_total_count=15,
    large_n=10,
    min_prop=0.7,
)
dc.pp.filter_by_prop(
    adata=pdata,
    min_prop=0.1,
    min_smpls=2,
)
```


```python
pdata.write_h5ad(snakemake.output.pseudo_rna, compression='gzip')
```

### ATAC pseudobulking




```python
# Open atac
adata = sc.read_h5ad(snakemake.input.merged_atac_anndata)

# Get pseudo-bulk profile, summing all raw counts by cell type and sample
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=snakemake.params.sample_key,
    groups_col=snakemake.params.separating_cluster,
    mode='sum'
)

# Filter out samples that have 
dc.pp.filter_samples(pdata, min_cells=snakemake.params.min_cells)

# Store raw counts in layers
pdata.layers["counts"] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata)
sc.tl.pca(pdata)

# Return raw counts to X
dc.pp.swap_layer(adata=pdata, key="counts", inplace=True)

# Save pseudobulked atac data
pdata.write_h5ad(snakemake.output.pseudo_atac, compression='gzip')
```

### chromVAR pseudobulking




```python
# Read in MuData
mdata = mu.read_h5mu(snakemake.input.merged_multiome)

# Separate out the chromvar modality
chromvar  = mdata.mod['chromvar']

# Add possibly missing data
for parameter in [snakemake.params.sample_key, 'celltype']:
    paramdict = mdata.obs[parameter].to_dict()
    chromvar.obs[parameter] = [paramdict[x] for x in chromvar.obs_names]

chromvar.obs['Sample-celltype'] = chromvar.obs[snakemake.params.sample_key].astype(str) + '_' + chromvar.obs[snakemake.params.separating_cluster].astype(str)
pchromvar = sc.get.aggregate(
    chromvar, 
    by='Sample-celltype', 
    func = 'mean'
    )
pchromvar.X = pchromvar.layers['mean']

pchromvar.write_h5ad(snakemake.output.pseudobulk_chromvar, compression='gzip')
```

## Differential analysis

### Differentially expressed genes




```python
# Import pseudobulked data as an AnnData object
# This needs to be Sample-Celltype by gene
pb = sc.read_h5ad(snakemake.input.pseudo_rna)

# Parameters to slice 
celltype_param = snakemake.params.celltype_params
cell_grouping = snakemake.params.celltypes
diagnosis_param = snakemake.params.diagnosis_param
diagnosis_control = snakemake.params.diagnosis_control
celltype = snakemake.params.celltype

# Change parameter names from the Decoupler hard-coded values to
# Dreampy hard-coded values
pb.obs['n_cells'] = pb.obs['psbulk_cells']
pb.obs['assays'] = pb.obs[celltype_param]

# Preprocessing
pb = dp.filter_samples(pb, min_cells=10, min_samples=3)
pb = dp.compute_tmm_factors(pb, assay_col=celltype_param)
assays_dict = dp.filter_by_expr(pb, assay_col=celltype_param)
```


```python
comparison_assay_pb = assay_pb[assay_pb.obs[diagnosis_param].isin(comparison)]

# log2CPM
comparison_assay_pb = dp.log2cpm(comparison_assay_pb)

if comparison_assay_pb.shape[0] != 0:

    # Voom weights
    comparison_assay_pb = dp.estimate_weights(
        comparison_assay_pb, formula=snakemake.params.formula, n_jobs=-1,
    )

    # Mixed-model fitting
    fit = dp.fit_models(
        comparison_assay_pb,
        formula=snakemake.params.formula,
        reml=True,
        n_jobs=snakemake.threads,
        assay_name=celltype,
    )
    print(f"  fit: {fit}")

    # eBayes
    fit_eb = dp.ebayes(fit)
    print(f"  eBayes: s2_prior={fit_eb.s2_prior:.4f}, "
        f"df_prior={fit_eb.df_prior:.1f}")

    # Results table
    results = dp.get_results(fit_eb, assay_name=celltype)
    results[celltype_param] = celltype
    results[diagnosis_param] = disease_name

    all_results = pd.concat([all_results, results])
```

### Differentially enriched motifs




```python
# Load RNA data
rna_adata_df = pd.read_csv('../../data/rna_cell_annot.csv')
rna_adata_df.index = rna_adata_df['Unnamed: 0']
del rna_adata_df['Unnamed: 0']
rna_adata_df

chromvar  = mu.read(snakemake.input.merged_multiome + "/mod/chromvar")


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
```

### Gene-set enrichment




```python
adata = sc.read_h5ad(snakemake.input.adata_path)

adata = adata[adata.obs[snakemake.params.separating_cluster] == snakemake.params.cell_type]

# Filter based on condition to test
diseaes_adata = adata[adata.obs[snakemake.params.disease_param].isin([snakemake.params.control, snakemake.params.disease])]

# Add one-hot column for disease comparison, sort the AnnData object
diseaes_adata.obs['disease'] = pd.Categorical(
    diseaes_adata.obs[snakemake.params.disease_param], 
    categories=[snakemake.params.disease, snakemake.params.control], 
    ordered=True)
indices = diseaes_adata.obs.sort_values(['disease']).index
diseaes_adata = diseaes_adata[indices,:].copy()
```


```python
# Import same ontology dataset as used in GREAT, reformat for GSEA
greatpy_go_set_df = pd.read_csv(snakemake.input.ontologies, delimiter=';')
greatpy_go_set_df['key'] = greatpy_go_set_df['name'] + ' (' + greatpy_go_set_df['id'] + ')'
go_mf = greatpy_go_set_df[['key', 'symbol']].groupby('key')['symbol'].apply(list).to_dict()

# Calculate GSEA enrichment
res = gp.gsea(
    data=diseaes_adata.to_df().T, # row -> genes, column-> samples
    gene_sets = go_mf,
    classes=[snakemake.params.disease, snakemake.params.control],
    permutation_num=1000,
    permutation_type='phenotype',
    outdir=None,
    method='log2_ratio_of_classes', # signal_to_noise
    threads= snakemake.threads # Speed up enrichment analysis
    )
```


```python
res.res2d
```

### Genomic Regions Enrichment of Annotations Tool (GREAT)




```python
# Create regulatory domain
regdom_hg38 = great.tl.create_regdom(
    tss_file = snakemake.input.tss_file,
    chr_sizes_file = snakemake.input.chr_sizes_file,
    association_rule = "basal_plus_extention",
    out_path = None
    )

# Load in differentially accessible region file
cell_type_DARs = pd.read_csv(snakemake.input.DAR_path)
```


```python
# Return list of top and bottom 5% enriched peaks
DAR_counts = cell_type_DARs.shape[0] # Total number of peaks
top_05_num = int(DAR_counts * 0.05) # 5% of the total number

cell_sign_DAR_df = cell_type_DARs.sort_values('log2FoldChange', ascending=False).iloc[:top_05_num]

# Split out chromosome, start and end for GREAT
cell_sign_DAR_df['chr'] = [x.split(':')[0] for x in cell_sign_DAR_df['peak']]
cell_sign_DAR_df['start'] = [x.split(':')[1].split('-')[0] for x in cell_sign_DAR_df['peak']]
cell_sign_DAR_df['end'] = [x.split(':')[1].split('-')[1] for x in cell_sign_DAR_df['peak']]

# Save split DataFrame
cell_sign_DAR_df[['chr', 'start', 'end']].to_csv(
    snakemake.output.cell_disease_peaks,
    sep = '\t',
    header=False,
    index=False
    )

# Compute GREAT enrichment
enrichment = great.tl.enrichment(
    test_file = snakemake.output.cell_disease_peaks,
    regdom_file = regdom_hg38,
    chr_size_file = snakemake.input.chr_sizes_file,
    annotation_file = snakemake.input.annotation_file,
    binom = True,
    hypergeom = True,
    )
```


```python

# Calculate Bonferroni FDR
enrichment = great.tl.set_fdr(enrichment, alpha = 0.01)
enrichment = great.tl.set_bonferroni(enrichment, alpha = 0.01)

# Save GREAT-based pathway enrichment
enrichment.to_csv(snakemake.output.cell_disease_GREAT, index=False)
```

## Correlation analysis

### Gene peak linkage




```python
conditions = snakemake.params.conditions
celltype = snakemake.params.celltype

# Gene names
gene_info_df = pd.read_csv(
    snakemake.input.gene_info, 
    delimiter='\t',
    names = ['ensembl', 'gene name', 'gene type'])

# TSS from Ensembl values
ensembl_TSS_df = pd.read_csv(
    snakemake.input.gene_tss,
    delimiter='\t',
    header=None,
    names=['Chr', 'Start', 'End', 'ensembl', '', 'strand'])
# Add the gene name for each start site
gene_TSS_df = pd.merge(
    left = ensembl_TSS_df[['Chr', 'Start', 'End', 'ensembl', 'strand']],
    right = gene_info_df[['ensembl', 'gene name']]
)
# Chromosome set
chromosomes = ['chr'+str(i+1) for i in range(22)] + ['chrX', 'chrY']

```


```python
# Load pseudobulked dataset
pdata = sc.read_h5ad(snakemake.input.pseudobulked_rna)
celltype_pdata = pdata[pdata.obs['celltype'] == celltype]
# Observation index needs to be stripped of 'SampleID_celltype'
celltype_pdata.obs_names = ['_'.join(x.split('_')[:-1]) for x in celltype_pdata.obs_names]

print('Getting list of co-accessible peaks')
# Dictionary to store pseudobulked chromatin accessibility datasets
celltype_condition_atac_dict = {}
```


```python
# Read in atac anndata object for each sample, 
for condition, temp_atac in zip(snakemake.params.conditions, snakemake.input.atac_files):
    print(f'Working on {condition} atac file')
    pdata_atac = dc.pp.pseudobulk(
        sc.read_h5ad(temp_atac),
        sample_col='sample_id',
        groups_col='sample_id',
        mode='sum'
    )
    # Store sum of counts
    pdata_atac.layers['counts'] = pdata_atac.X.copy()

    # Log-normalized
    # Normalize and scale
    sc.pp.normalize_total(pdata_atac, target_sum=1e4)
    del temp_atac
    celltype_condition_atac_dict[condition] = pdata_atac

```


```python
print('Get overlapping peaks')

# Define file names
file_names = snakemake.input.bed_files
# Read in files
celltype_beds = [pr.read_bed(x) for x in file_names]
# Concatenate and merge all overlapping peaks
overlapping_peaks = pr.concat(celltype_beds).merge()
# Add columns with overlapping boundaries
peakset_overlap_df = pr.concat([x.join(overlapping_peaks, suffix = '_overlap') for x in celltype_beds]).df
peakset_overlap_df['overlap peak'] = peakset_overlap_df['Chromosome'].astype(str) + ':' + peakset_overlap_df['Start_overlap'].astype(str) + '-' + peakset_overlap_df['End_overlap'].astype(str)
peakset_overlap_df['peak'] = peakset_overlap_df['Chromosome'].astype(str) + ':' + peakset_overlap_df['Start'].astype(str) + '-' + peakset_overlap_df['End'].astype(str)
peakset_overlap = peakset_overlap_df[['peak','overlap peak']]
peak2overlap = dict(zip(peakset_overlap['peak'], peakset_overlap['overlap peak']))

```


```python
print('Creating dictionary of co-accessible peaks')
# Dictionary to store results in 
celltype_condition_coacc_dict = {}
for condition, circe_file in zip(conditions, snakemake.input.circe_files):
    control_network = pd.read_csv(
        circe_file,
        delimiter='\t', 
        header=None, 
        names = ['chr 1', 'start 1', 'end 1', 'chr 2', 'start 2', 'end 2', 'score'])
    control_network = control_network[control_network['score'] > 0.2]
    celltype_condition_coacc_dict[condition] = control_network

```


```python
print('Create promoter-enhancer list for each condition')
# DataFrame to store values
TSS_coaccess_df = pd.DataFrame()
# Iterate over condtions
for condition in conditions:
    print(f'Working on condition {condition}')
    # Select network
    circe_network = celltype_condition_coacc_dict[condition]

    # Iterate over chromosomes
    for chrom in chromosomes:
        print(chrom)
        test_cross = pd.merge(
            left = circe_network[circe_network['chr 1'] == chrom].sort_values('start 1'),
            right = gene_TSS_df[gene_TSS_df['Chr'] == chrom],
            how='cross'
            )
        promoter_coaccessible_first_peak = test_cross[
            ((test_cross['start 1'] < test_cross['Start']-200) & (test_cross['end 1'] > test_cross['Start']-200)) |
            ((test_cross['start 1'] < test_cross['Start']+200) & (test_cross['end 1'] > test_cross['Start']+200))
        ]
        promoter_coaccessible_first_peak['promoter'] = 'peak 1'

        promoter_coaccessible_second_peak = test_cross[
            ((test_cross['start 2'] < test_cross['Start']-200) & (test_cross['end 2'] > test_cross['Start']-200)) |
            ((test_cross['start 2'] < test_cross['Start']+200) & (test_cross['end 2'] > test_cross['Start']+200))
        ]
        promoter_coaccessible_second_peak['promoter'] = 'peak 2'
        promoter_coaccessible = pd.concat([promoter_coaccessible_first_peak, promoter_coaccessible_second_peak])

        promoter_coaccessible['diagnosis'] = condition
        TSS_coaccess_df = pd.concat([TSS_coaccess_df, promoter_coaccessible])
TSS_coaccess_df['peak 1'] = TSS_coaccess_df['chr 1'] + ':' + TSS_coaccess_df['start 1'].astype(str) + '-' + TSS_coaccess_df['end 1'].astype(str)
TSS_coaccess_df['peak 2'] = TSS_coaccess_df['chr 2'] + ':' + TSS_coaccess_df['start 2'].astype(str) + '-' + TSS_coaccess_df['end 2'].astype(str)
TSS_coaccess_df = TSS_coaccess_df[['peak 1', 'peak 2', 'Start', 'score', 'gene name', 'diagnosis', 'promoter']]

```


```python
# Find overlap peak key when enhancer is peak 2
peak1enhancer_df = TSS_coaccess_df[TSS_coaccess_df['promoter'] == 'peak 2']
peak1enhancer_df['overlap promoter peak'] = [peak2overlap[x] for x in peak1enhancer_df['peak 2']]
peak1enhancer_df['overlap enhancer peak'] = [peak2overlap[x] for x in peak1enhancer_df['peak 1']]

# Find overlap peak key when enhancer is peak 2
peak2enhancer_df = TSS_coaccess_df[TSS_coaccess_df['promoter'] == 'peak 1']
peak2enhancer_df['overlap promoter peak'] = [peak2overlap[x] for x in peak2enhancer_df['peak 1']]
peak2enhancer_df['overlap enhancer peak'] = [peak2overlap[x] for x in peak2enhancer_df['peak 2']]

# Add overlap peak for enhancer
TSS_coaccess_df = pd.concat([peak1enhancer_df, peak2enhancer_df])

```


```python
print('Promoter enhancer list is done')

# Dictionary of pseudobulks merged
rna_atac_merge_dict = {x: pd.merge(left = celltype_condition_atac_dict[x].to_df(),right = celltype_pdata.to_df(), left_index = True, right_index = True) for x in ['control', 'PD', 'LBD']}

p_vals= []
stats = []
print('Computing stats for each promoter-enhancer')
for i in range(TSS_coaccess_df.shape[0]):
    temp_df = TSS_coaccess_df.iloc[[i]]
    if temp_df['promoter'].values[0] == 'peak 1':
        peak_name = 'peak 2'
    else:
        peak_name = 'peak 1'
    temp_condition = temp_df['diagnosis'].values[0]
    peak = temp_df[peak_name].values[0]
    gene = temp_df['gene name'].values[0]
    if gene in celltype_pdata.var_names:
        stat, p = scipy.stats.spearmanr(rna_atac_merge_dict[temp_condition][gene], rna_atac_merge_dict[temp_condition][peak])
    else:
        stat, p = 0, 1
    p_vals.append(p)
    stats.append(stat)
TSS_coaccess_df['gene-peak link stat'] = stats
TSS_coaccess_df['gene-peak link p-value'] = p_vals

TSS_coaccess_df.to_csv(snakemake.output.gene_peak_linkage, index=False)
```

### Gene motif linkage




```python
# Import the celltype of interest
celltype_param = snakemake.params.celltype_param
celltype = snakemake.params.celltype
# Import the diagnosis of interest, and the parameter that the DataFrame is sliced from
diagnosis_param = snakemake.params.diagnosis_param
diagnosis = snakemake.params.diagnosis
```


```python
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
```


```python
# Save the motif names to a list, convert them be separated by _
TF_motif_names = pchromvar_df.columns.to_list()
TF_motif_names = ['_'.join(x.replace("::", "..").replace("-", ".").split('.')) for x in TF_motif_names]
# Rename the columns with the corrected TF motif signals
pchromvar_df = pchromvar_df.rename(columns=dict(zip(pchromvar_df.columns.to_list(), TF_motif_names)))

# Combine the RNA and TF motif signals
prna_df = prna.to_df()
gene_names = prna_df.columns.to_list()
```


```python
# Add new scripts
gene_motif_df = pd.merge(
    left = prna_df,
    right = pchromvar_df,
    left_index = True,
    right_index = True
)
gene_motif_df = pd.merge(
    left = gene_motif_df,
    right = prna.obs[[celltype_param, diagnosis_param]],
    left_index = True,
    right_index = True
)
```


```python
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
diagnosis_gene_motif_df['r-squared'] = diagnosis_gene_motif_df['Pearson r'] ** 2
```


```python
# Export gene-motif linkages
diagnosis_gene_motif_df.to_csv(snakemake.output.gene_motif_links, index = False)
```

### Cell cell communication




```python
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

control_adata = adata[adata.obs[snakemake.params.disease_param] == snakemake.params.control]

# Run rank_aggregate
li.mt.cellphonedb(
    control_adata,
    groupby='celltype',
    #resource=cellphone_resource,
    resource_name='consensus',
    expr_prop=0.1,
    verbose=True,
    key_added='li_res',
    n_jobs=snakemake.threads
    )

li.mt.rank_aggregate(control_adata, 
                     groupby='celltype',
                     resource_name='consensus',
                     expr_prop=0.1,
                     verbose=True)

cellphone_df = control_adata.uns['li_res'].copy()
cellphone_df.to_csv(snakemake.output.cell_cell_communication_data, index=None)
```

### Co-accessibility




```python
adata = sc.read_h5ad(snakemake.input.celltype_atac)

ci.add_region_infos(adata, sep=(':', '-'))
adata.var['start'] = adata.var['start'].astype(int)
adata.var['end'] = adata.var['end'].astype(int)
adata.var['chr'] = adata.var['chromosome']

# Convert
adata = ci.metacells.compute_metacells(
    adata, 
    method='sum', 
    k=25, 
)
```


```python
# Get co-accessibility scores
final_score = ci.sliding_graphical_lasso(
    adata,
    n_samples=50,
    n_samples_maxtry=100,
    max_alpha_iteration=500,
    verbose=True,
    njobs=snakemake.threads
)

adata.varp['atac_network'] = final_score

# Compute the co-accessibility network
ci.compute_atac_network(adata, njobs=snakemake.threads)

adata.write_h5ad(snakemake.output.celltype_atac, compression='gzip')
```


```python
# Extract the network and find CCANs modules
circe_network = ci.extract_atac_links(adata)

circe_network['Peak1'] = circe_network['Peak1'].str.replace(':', '_').str.replace('-', '_')
circe_network['Peak2'] = circe_network['Peak2'].str.replace(':', '_').str.replace('-', '_')

# Split columns for links files
circe_network['chr 1'] = [x.split('_')[0] for x in circe_network['Peak1']]
circe_network['start 1'] = [x.split('_')[1] for x in circe_network['Peak1']]
circe_network['end 1'] = [x.split('_')[2] for x in circe_network['Peak1']]

circe_network['chr 2'] = [x.split('_')[0] for x in circe_network['Peak2']]
circe_network['start 2'] = [x.split('_')[1] for x in circe_network['Peak2']]
circe_network['end 2'] = [x.split('_')[2] for x in circe_network['Peak2']]

links_network = circe_network[['chr 1', 'start 1', 'end 1', 'chr 2', 'start 2', 'end 2', 'score']]
links_network.to_csv(snakemake.output.circe_network, sep='\t', index=None, header=None)
```
