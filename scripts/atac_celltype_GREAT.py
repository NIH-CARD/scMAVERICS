import numpy as np
import pandas as pd

import greatpy as great

# Import celltype
celltype = snakemake.params.celltype

# Import consensus bed file
consensus_peaks = pd.read_csv(
    snakemake.input.consensus_bed,
    delimiter = '\t',
    header = None,
    names = ['Chr', 'Start', 'End', 'source peaks', 'score', '.']
    )
consensus_peaks['number of source peaks'] = [len(x.split(',')) for x in consensus_peaks['source peaks']]
consensus_peaks['peak'] = consensus_peaks['Chr'].astype(str) + ':' + consensus_peaks['Start'].astype(str) + '-' + consensus_peaks['End'].astype(str)

for celltype in snakemake.params.cell_types:
    consensus_peaks[celltype] = consensus_peaks['source peaks'].str.contains(celltype)

def celltype_peak_filter(included_celltypes, excluded_celltypes=[]):
    new_consensus_peaks = consensus_peaks

    for celltype in included_celltypes:
        new_consensus_peaks = new_consensus_peaks[new_consensus_peaks[celltype]]
    if len(excluded_celltypes) > 0:
        for celltype in excluded_celltypes:
            new_consensus_peaks = new_consensus_peaks[~(new_consensus_peaks[celltype])]
    return new_consensus_peaks

missing_one_celltype = [x for x in cell_types if x != celltype]
celltype_peaks = celltype_peak_filter(included_celltypes=[celltype], excluded_celltypes=missing_one_celltype)

# Create regulatory domain
regdom_hg38 = great.tl.create_regdom(
    tss_file = snakemake.input.tss_file,
    chr_sizes_file = snakemake.input.chr_sizes_file,
    association_rule = "basal_plus_extention",
    out_path = None
    )

# Compute GREAT enrichment
enrichment = great.tl.enrichment(
    test_file = celltype_peaks,
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
