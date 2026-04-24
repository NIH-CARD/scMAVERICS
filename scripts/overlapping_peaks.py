import pyranges as pr
import pandas as pd

# Define file names
files_names = snakemake.input.peak_files
#file_names = [f'../../data/celltypes/{celltype}/{celltype}_{condition}_peaks.bed' for condition in conditions]
# Read in files
celltype_beds = [pr.read_bed(x) for x in file_names]
# Concatenate and merge all overlapping peaks
overlapping_peaks = pr.concat(celltype_beds).merge()

# Save overlapping peaks as a bed file
overlapping_peaks.to_bed(snakemake.output.celltype_overlapping_peaks)

# Add columns with overlapping boundaries
peakset_overlap_df = pr.concat([x.join(overlapping_peaks, suffix = '_overlap') for x in celltype_beds]).df
peakset_overlap_df['overlap peak'] = peakset_overlap_df['Chromosome'].astype(str) + ':' + peakset_overlap_df['Start_overlap'].astype(str) + '-' + peakset_overlap_df['End_overlap'].astype(str)
peakset_overlap_df['peak'] = peakset_overlap_df['Chromosome'].astype(str) + ':' + peakset_overlap_df['Start'].astype(str) + '-' + peakset_overlap_df['End'].astype(str)

# Save all data in csv
peakset_overlap_df.to_csv(snakemake.output.celltype_overlapping_celltype_peaks, index=False)
