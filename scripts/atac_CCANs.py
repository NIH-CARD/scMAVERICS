import numpy as np
import pandas as pd
import scanpy as sc
import circe as ci

# Load in cell type
cell_type_atac = sc.read_h5ad(snakemake.input.celltype_atac)

# Get previously calculated network
circe_network = ci.extract_atac_links(cell_type_atac)

# Calculate CCAN
ccans = ci.find_ccans(circe_network, seed=0)

# Assign each peak with a CCAN, if it is part of a network
ccans['Chr'] = [x.split(':')[0] for x in ccans['Peak']]
ccans['Start'] = [int(x.split(':')[1].split('-')[0]) for x in ccans['Peak']]
ccans['End'] = [int(x.split(':')[1].split('-')[1]) for x in ccans['Peak']]

ccans.to_csv(snakemake.output.output_DAR_CCAN_data , index=False)