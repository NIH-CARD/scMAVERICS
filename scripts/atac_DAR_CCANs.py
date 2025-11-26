import numpy as np
import pandas as pd
import scanpy as sc
import circe as ci

# Load in cell type
cell_type_atac = sc.read_h5ad(snakemake.input.celltype_atac)

# Get DARs 
cell_type_DARs = pd.read_csv(snakemake.input.output_DAR_data)

# Get previously calculated network
circe_network = ci.extract_atac_links(cell_type_atac)

# Calculate CCAN
ccans = ci.find_ccans(circe_network, seed=0)

# For cell types with 
peak2ccan = dict(zip(ccans['Peak'], ccans['CCAN']))

# Assign each peak with a CCAN, if it is part of a network
ccan_assign = []
for x in cell_type_DARs['peak']:
    if x in peak2ccan.keys():
        ccan_assign.append(peak2ccan[x])
    else:
        ccan_assign.append(None)

cell_type_DARs['CCAN'] = ccan_assign

cell_type_DARs.to_csv(snakemake.output.output_DAR_CCAN_data , index=False)