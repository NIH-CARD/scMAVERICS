import pandas as pd
import pyranges as pr

# Load bed file
peak_file = pd.read_csv(snakemake.input.xls, delimiter='\t', skiprows=21)

# Load in blacklist
black_list = pr.read_bed(snakmake.input.blacklist)

# Convert to bed file
peak_pr = pr.PyRanges(
    chromosomes=peak_file['chr'],
    starts=peak_file['start'],
    ends=peak_file['end']
    )

acceptable_peaks = ['chr' + str(i+1) for i in range(0, 22)] + ['chrX', 'chrY']

filter_intermed = peak_pr.join(black_list, how='outer')
filter_intermed = filter_intermed[filter_intermed.df['Start_b'] == -1]
peak_pr = filter_intermed[filter_intermed.df['Chromosome'].isin(acceptable_peaks)]


# Save the bed file
peak_pr[['Chromosome', 'Start', 'End']].to_bed(snakemake.output.cell_bedfile)
