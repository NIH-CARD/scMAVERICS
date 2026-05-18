import pandas as pd
import scanpy as sc
import numpy as np
import scipy

""" Define analysis parameters """
celltype = snakemake.params.cell_type
diagnosis = snakemake.params.disease
sample_param = snakemake.params.sample_key
disease_param = snakemake.params.disease_param

""" Load in genome-specific TSS sites """

# Gene names
gene_info_df = pd.read_csv(
    snakemake.input.gene_into, 
    delimiter='\t',
    names = ['ensembl', 'gene name', 'gene type'])

# TSS from Ensembl values
ensembl_TSS_df = pd.read_csv(
    snakemake.input.tss_file,
    delimiter='\t',
    header=None,
    names=['Chr', 'Start', 'End', 'ensembl', '', 'strand'])
# Add the gene name for each start site
gene_TSS_df = pd.merge(
    left = ensembl_TSS_df[['Chr', 'Start', 'End', 'ensembl', 'strand']],
    right = gene_info_df[['ensembl', 'gene name']]
)

# Chromosome set (human)
chromosomes = ['chr'+str(i+1) for i in range(22)] + ['chrX', 'chrY']

""" Load in pseudobulked RNA values and add conditions"""

# Load pseudobulked dataset 
pdata = sc.read_h5ad(snakemake.input.pseudo_rna)

# Create DataFrame to quickly pull correlation values
pdata_df = pdata.to_df()
pdata_df['celltype'] = [x.split('_')[-1] for x in pdata_df.index]
pdata_df['Sample_ID'] = ['_'.join(x.split('_')[:-1]) for x in pdata_df.index]

# Convert diagnosis parameter 
sample_diagnosis = dict(zip(pdata.obs['Sample_ID'], pdata.obs[disease_param]))
pdata_df['diagnosis'] = [sample_diagnosis[x] for x in pdata_df['Sample_ID']]

# Import CCAN values
ccan_df = pd.read_csv(snakemake.input.output_CCAN_data)

""" Find all genes whose promoter lies within a CCAN """

# Initialize DataFrame to hold which genes are in which CCAN, chromosomes
ccan_gene_df = pd.DataFrame()

# Iterate over chromosomes (faster TSS filtering using the pd.merge with how='cross')
for chrom in chromosomes:

    # Filter for condition
    temp_filt_CCAN_df =  ccan_df[ccan_df['Chr'] == chrom]

    # Overlap all rows (easier search)
    test_cross = pd.merge(
        left = temp_filt_CCAN_df,
        right = gene_TSS_df[gene_TSS_df['Chr'] == chrom],
        how='cross',
        suffixes = ['', '_TSS']
        )

    # Filter DataFrame for TSS within +/- 200bp of peak
    overlapping_peaks = test_cross[
        ((test_cross['Start'] < test_cross['Start_TSS']-200) & (test_cross['End'] > test_cross['Start_TSS']-200)) |
        ((test_cross['Start'] < test_cross['Start_TSS']+200) & (test_cross['End'] > test_cross['Start_TSS']+200))
    ]

    # Assign cell type and diagnosis values
    overlapping_peaks['celltype'] = celltype
    overlapping_peaks['diagnosis'] = diagnosis

    # Add to the overall CCAN_df
    ccan_gene_df = pd.concat([ccan_gene_df, overlapping_peaks[['celltype', 'diagnosis', 'Chr', 'CCAN', 'gene name']].drop_duplicates()])

# Save the CCAN gene hub
ccan_gene_df.to_csv(snakemake.output.ccan_gene, index=False)

""" Co-expression of genes from intra-CCANs """
# Get list of all genes in this condition
condition_pdata_df = pdata_df[(pdata_df['celltype'] == celltype) & (pdata_df['diagnosis'] == diagnosis)]

# Slice out data
ccan_condition = ccan_gene_df[
    (ccan_gene_df['celltype'] == celltype) & 
    (ccan_gene_df['diagnosis'] == diagnosis)]

# Initialize CCAN data list
correlation_data = []

# List of CCANs to iterage through
ccan_num = ccan_condition.CCAN.drop_duplicates().to_list()
for ccan in ccan_num:

    # Slice out that CCAN
    specific_ccan = ccan_condition[ccan_condition['CCAN'] == ccan]

    # Iterate through all genes in that CCAN, making sure not to count twice
    for i, gene_1 in enumerate(specific_ccan['gene name']):
        for j, gene_2 in enumerate(specific_ccan['gene name']):
            if i < j and gene_1 in condition_pdata_df.columns and gene_2 in condition_pdata_df.columns:

                # Compute Spearman correlation
                rho, p = scipy.stats.spearmanr(condition_pdata_df[gene_1], condition_pdata_df[gene_2])
                correlation_data.append([celltype, diagnosis, gene_1, gene_2, rho, p])

# Convert data to DataFrame
correlation_data_df = pd.DataFrame(
    correlation_data, 
    columns=['celltype', 'diagnosis', 'gene 1', 'gene 2', 'Spearman statistic', 'Spearman p-value']
    )

""" Co-expression of genes from inter-CCANs """

# Initialize list to store data 
non_ccan_correlation_data = []

# Only sample from CCANs that have more than one gene
ccan_count_df = ccan_condition.groupby(['Chr', 'CCAN']).count().reset_index()
good_CCANs = ccan_count_df[ccan_count_df['celltype'] > 2].CCAN.values

# Only sample from chromosomes that have more than one CCAN
chr_count_df = ccan_count_df[ccan_count_df['celltype'] > 2].groupby('Chr').count().reset_index()
good_Chrs = chr_count_df[chr_count_df['celltype'] > 1].Chr.values

# Filter out data 
ccan_condition = ccan_condition[
    (ccan_condition['gene name'].isin(condition_pdata_df.columns)) &
    (ccan_condition['CCAN'].isin(good_CCANs)) &
    (ccan_condition['Chr'].isin(good_Chrs))
    ]

# For reproducible randomness, import seed
np.random.seed(snakemake.params.random_seed)

# List of all genes to sample from 
all_genes = ccan_condition['gene name'].to_list()

# Initialize looping counter, to make sure same number of random inter-comparisons as intra-comparisons
k = 0
while k < len(correlation_data_df):

    # Randomly pick first gene
    gene_1 = np.random.choice(all_genes)

    # Randomly pick second gene that is in the same chromosome and not in the same CCAN
    gene_2 = np.random.choice(ccan_condition[
        (ccan_condition['Chr'] == ccan_condition[ccan_condition['gene name'] == gene_1].Chr.values[0]) &
        (ccan_condition['CCAN'] != ccan_condition[ccan_condition['gene name'] == gene_1].CCAN.values[0])]['gene name'])
    
    # Calculate the Spearman statistic
    rho, p = scipy.stats.spearmanr(condition_pdata_df[gene_1], condition_pdata_df[gene_2])
    non_ccan_correlation_data.append([celltype, diagnosis, gene_1, gene_2, rho, p])
    
    # Add to counter
    k += 1

# Convert to DataFrame
non_correlation_data_df = pd.DataFrame(non_ccan_correlation_data, columns=['celltype', 'diagnosis', 'gene 1', 'gene 2', 'Spearman statistic', 'Spearman p-value'])

# Combined inter- and intra- correlation DataFrames, saving the correlation type
correlation_data_df['similarity'] = 'intra'
non_correlation_data_df['similarity'] = 'inter'
hub_coexpression_df = pd.concat([correlation_data_df, non_correlation_data_df])

# Export to .csv
hub_coexpression_df.to_csv(snakemake.output.ccan_corr, index=False)
