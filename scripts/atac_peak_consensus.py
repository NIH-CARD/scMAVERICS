import numpy as np
import pandas as pd
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.pseudobulk_peak_calling import peak_calling
import scanpy as sc
import pyranges
import os

# Read in narrow peak dictionaries
peak_files = snakemake.input.narrow_peaks

# List of cell types
cells = [x.split('/')[-1].split('_')[0] for x in peak_files]

# Create narrow peak object
narrow_peak_dict = dict(zip(cells, [pyranges.PyRanges(pyranges.read_bed(x).df.rename(columns={'ThickStart': 'FC_summit', 'ThickEnd': '-log10_pval', 'ItemRGB': '-log10_qval', 'BlockCount': 'Summit'}, copy=True)) for x in peak_files]))

# Read in chromosome lengths
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

# Set peak target width
peak_half_width=250

# Get consensus peaks from peak dictionary
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes)

# Export the consensus bed file
consensus_peaks.to_bed(
    path = snakemake.output.consensus_bed,
    keep =True,
    compression = 'infer',
    chain = False)