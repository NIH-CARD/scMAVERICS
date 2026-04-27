import numpy as np
import scanpy as sc
import scipy
import decoupler as dc
import pandas as pd
import pyranges as pr

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

# Load pseudobulked dataset
pdata = sc.read_h5ad(snakemake.input.pseudobulked_rna)
celltype_pdata = pdata[pdata.obs['celltype'] == celltype]
# Observation index needs to be stripped of 'SampleID_celltype'
celltype_pdata.obs_names = ['_'.join(x.split('_')[:-1]) for x in celltype_pdata.obs_names]

print('Getting list of co-accessible peaks')
# Dictionary to store pseudobulked chromatin accessibility datasets
celltype_condition_atac_dict = {}

# Read in atac anndata object for each sample, 
for condition, temp_atac in zip(conditions, atac_files):
    print(f'Working on {condition} atac file')
    pdata_atac = dc.pp.pseudobulk(
        temp_atac,
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

print('Get overlapping peaks')

# Define file names
file_names = [f'../../data/celltypes/{celltype}/{celltype}_{condition}_peaks.bed' for condition in conditions]
# Read in files
celltype_beds = [pr.read_bed(x) for x in file_names]
# Concatenate and merge all overlapping peaks
overlapping_peaks = pr.concat(celltype_beds).merge()
# Add columns with overlapping boundaries
peakset_overlap_df = pr.concat([x.join(overlapping_peaks, suffix = '_overlap') for x in celltype_beds]).df
peakset_overlap_df['overlap peak'] = peakset_overlap_df['Chromosome'].astype(str) + ':' + peakset_overlap_df['Start_overlap'].astype(str) + '-' + peakset_overlap_df['End_overlap'].astype(str)
peakset_overlap_df['peak'] = peakset_overlap_df['Chromosome'].astype(str) + ':' + peakset_overlap_df['Start'].astype(str) + '-' + peakset_overlap_df['End'].astype(str)

print('Creating dictionary of co-accessible peaks')
# Dictionary to store results in 
celltype_condition_coacc_dict = {}
for condition, circe_file in zip(conditions, circe_files):
    control_network = pd.read_csv(
        circe_file,
        delimiter='\t', 
        header=None, 
        names = ['chr 1', 'start 1', 'end 1', 'chr 2', 'start 2', 'end 2', 'score'])
    control_network = control_network[control_network['score'] > 0.2]
    celltype_condition_coacc_dict[condition] = control_network

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
TSS_coaccess_df = TSS_coaccess_df[['peak 1', 'peak 2', 'Start', 'gene name', 'diagnosis', 'promoter']]

# Find overlap peak key when enhancer is peak 2
peak1enhancer_df = TSS_coaccess_df[TSS_coaccess_df['promoter'] == 'peak 2']
peak1enhancer_df['overlap enhancer peak'] = [peak2overlap[x] for x in peak1enhancer_df['peak 1']]

# Find overlap peak key when enhancer is peak 2
peak2enhancer_df = TSS_coaccess_df[TSS_coaccess_df['promoter'] == 'peak 1']
peak2enhancer_df['overlap enhancer peak'] = [peak2overlap[x] for x in peak2enhancer_df['peak 2']]

# Add overlap peak for enhancer
TSS_coaccess_df = pd.concat([peak1enhancer_df, peak2enhancer_df])

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