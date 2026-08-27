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

# Filter for only the main 23 chromosomes
pseudo_fragment = pseudo_fragment.filter(pseudo_fragment['chrom'].is_in(pr.data.chromsizes().df['Chromosome'].to_list()))
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
    score=pseudo_fragment['score'].to_pandas()
)

chrom_sizes_dict = dict(zip(chromosizes['Chromosome'], chromosizes['length']))

bw = pyBigWig.open(snakemake.output.celltype_bigwig, "w")
bw.addHeader([(chrom, int(size)) for chrom, size in chrom_sizes_dict.items()])

df_sorted = combined_bed.df.sort_values(by=["Chromosome", "Start"])

# Stream entries chromosome by chromosome straight to disk
for chrom, group in df_sorted.groupby("Chromosome"):
    if chrom not in chrom_sizes_dict:
        continue
        
    bw.addEntries(
        [str(chrom)] * len(group),
        group["Start"].tolist(),
        ends=group["End"].tolist(),
        values=group["score"].astype(float).tolist() 
    )

combined_bed.to_bigwig(
    snakemake.output.celltype_normalized_bigwig,
    chromosome_sizes = dict(zip(chromosizes['Chromosome'], chromosizes['length'])),
    rpm=True)