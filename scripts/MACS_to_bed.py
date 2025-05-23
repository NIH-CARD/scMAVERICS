import pandas as pd
import pyranges as pr

# Load bed file
peak_file = pd.read_csv(snakemake.input.xls, delimiter='\t', skiprows=21)

# Convert to bed file
peak_pr = pr.PyRanges(
    chromosomes=peak_file['chr'],
    starts=peak_file['start'],
    ends=peak_file['end']
    )

# Save the bed file
peak_pr.to_bed(snakemake.output.cell_bedfile)
