# Import general packages
import numpy as np
import pandas as pd
import gc
# Import specialized packages
import scanpy as sc
import pychromvar as pc
import muon as mu
from pyjaspar import jaspardb
import scipy.sparse as sp

# Load MuData object
mdata = mu.read(snakemake.input.merged_multiome)

mdata['atac'].X = mdata['atac'].X.astype(np.float64)

# Get reference genome
pc.add_peak_seq(
    mdata['atac'],
    genome_file = snakemake.input.reference_genome,
    delimiter = ':|-'
)

# Correct for GC bias
global_bg_peaks = pc.add_gc_bias(mdata['atac'])
pc.get_bg_peaks(
    mdata['atac'],
    n_jobs = snakemake.threads
)

# Convert the legacy matrix to a modern sparse array format
mdata['atac'].X = sp.csr_array(mdata['atac'].X)
# Compute TF motif expected activity
global_expectations = pc.compute_expectation(mdata['atac'].X)

# Fetch motifs (modify if non-human)
jdb_obj = jaspardb(release='JASPAR2026')
motifs = jdb_obj.fetch_motifs(
    collection = 'CORE',
    species=9606)

# Match motifs in each sample
pc.match_motif(mdata['atac'], motifs=motifs)

# Return expectation values 
#dev = pc.compute_deviations(mdata, n_jobs = -1)

# Custom chunked method
chunk_size = 10000
total_cells = mdata['atac'].n_obs
chunks = [mdata['atac'][i:i + chunk_size].copy() for i in range(0, total_cells, chunk_size)]

dev = []

for i in range(0, mdata['atac'].n_obs, chunk_size):
    print(f"Processing chunk {i // chunk_size + 1}...")
    
    # Extract only the sparse chunk array to minimize memory footprint
    chunk_matrix = mdata['atac'].X[i:i + chunk_size, :].copy() 
    
    # Run sequentially (n_jobs=1) to prevent joblib from cloning memory
    dev_chunk = pc.compute_deviations(
        chunk_matrix, 
        bg_peaks=global_bg_peaks, 
        expectation=global_expectations,
        n_jobs=1
    )
    dev.append(dev_chunk)
    del chunk_matrix
    gc.collect()

# Save the motif value AnnData inside the muon object
mdata.mod['chromvar'] = dev
mdata.mod['chromvar'].raw = dev

# Write the muon object
mdata.write(snakemake.output.merged_multiome, compression='gzip')