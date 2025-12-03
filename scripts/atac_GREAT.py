import numpy as np
import pandas as pd

import greatpy as great

# Create regulatory domain
regdom_hg38 = great.tl.create_regdom(
    tss_file = snakemake.input.tss_file,
    chr_sizes_file = snakemake.input.chr_sizes_file,
    association_rule = "basal_plus_extention",
    out_path = None
    )

# Load in differentially accessible region file
cell_type_DARs = pd.read_csv(snakemake.input.DAR_path)

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

# Calculate Bonferroni FDR
enrichment = great.tl.set_fdr(enrichment, alpha = 0.01)
enrichment = great.tl.set_bonferroni(enrichment, alpha = 0.01)

# Save GREAT-based pathway enrichment
enrichment.to_csv(snakemake.output.cell_disease_GREAT, index=False)