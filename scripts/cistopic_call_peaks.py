import numpy as np
import pandas as pd
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.pseudobulk_peak_calling import peak_calling
import scanpy as sc
import pickle 
import os

# Load in the locations of the bed files, per cell type
bed_paths = {}
with open(snakemake.input.bed_paths) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

# Create directory to store peak-calling data
os.makedirs(snakemake.params.MACS_dir, exist_ok = True)

# Call peaks from bed paths
narrow_peak_dict = peak_calling(
    macs_path = "macs2",
    bed_paths = bed_paths,
    outdir = snakemake.params.MACS_dir,
    genome_size = 'hs',
    n_cpu = 32,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.01,
    nolambda = False,
    skip_empty_peaks=True,
    _temp_dir = snakemake.resources.tmpdir
)

# Save peaks as dictionary
with open(snakemake.output.peak_dict, 'wb') as f:
  pickle.dump(narrow_peak_dict, f)

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