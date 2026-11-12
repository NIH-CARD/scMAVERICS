import sys
import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr
import pyBigWig

pseudo_fragment = pl.read_csv(
    snakemake.input.pseudo_fragment_file, 
    separator='\t', 
    comment_prefix='#', 
    n_threads=snakemake.threads)
pseudo_fragment.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']

# Convert to pyranges
combined_bed = pr.PyRanges(
    chromosomes = pseudo_fragment['chrom'],
    starts=pseudo_fragment['chromStart'],
    ends=pseudo_fragment['chromEnd']
    )

# Read in Chromosome sizes
chromosizes = pd.read_csv(
    snakemake.input.genome_length, 
    header=None, 
    delimiter='\t',
    names=['Chromosome', 'length']
    )

combined_bed = pr.PyRanges(
    chromosomes=pseudo_fragment['chrom'].to_pandas(),
    starts=pseudo_fragment['chromStart'].to_pandas(),
    ends=pseudo_fragment['chromEnd'].to_pandas(),
)

chrom_sizes_dict = dict(zip(chromosizes['Chromosome'], chromosizes['length']))

combined_bed.to_bigwig(
    snakemake.output.celltype_bigwig,
    chromosome_sizes = chrom_sizes_dict,
    rpm=False)

combined_bed.to_bigwig(
    snakemake.output.celltype_normalized_bigwig,
    chromosome_sizes = chrom_sizes_dict,
    rpm=True)