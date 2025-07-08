import sys
import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr


pseudo_fragment = pl.read_csv(snakemake.input.pseudo_fragment_file, separator='\t', comment_prefix='#', n_threads=8)
pseudo_fragment.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']

# Filter for only the main 23 chromosomes
pseudo_fragment = pseudo_fragment.filter(pseudo_fragment['chrom'].is_in(pr.data.chromsizes().df['Chromosome'].to_list()))
# Convert to pyranges
combined_bed = pr.PyRanges(
    chromosomes = pseudo_fragment['chrom'],
    starts=pseudo_fragment['chromStart'],
    ends=pseudo_fragment['chromEnd']
    )
combined_bed.to_bigwig(snakemake.output.celltype_bigwig, rpm=False)
combined_bed.to_bigwig(snakemake.output.celltype_normalized_bigwig, rpm=True)