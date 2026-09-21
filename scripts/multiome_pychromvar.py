import numpy as np
import pandas as pd
import gc
import scanpy as sc
import muon as mu
import pychromvar as pc
from pyjaspar import jaspardb
import scipy.sparse as sp

# Load MuData object
mdata = mu.read(snakemake.input.merged_multiome)

# CRITICAL: Keep memory light using float32 instead of float64
mdata['atac'].X = mdata['atac'].X.astype(np.float32)

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

jdb_obj = jaspardb(release='JASPAR2026')
motifs = jdb_obj.fetch_motifs(collection='CORE', species=9606)

pc.match_motif(mdata['atac'], motifs=motifs)

# Chunk 
chunk_size = snakemake.params.chunk_size  # Working in chunks of 100000; adjust based on RAM
total_cells = mdata['atac'].n_obs
dev_chunks = []

for i in range(0, total_cells, chunk_size):
    print(f"Processing chunk {i // chunk_size + 1}...")
    
    # Slice the AnnData object in-place without cloning the background data
    chunk_adata = mdata['atac'][i:i + chunk_size].copy()
    
    # compute_deviations usually expects an AnnData object containing the motif match annotations
    dev_chunk = pc.compute_deviations(
        chunk_adata,
        n_jobs = 1  # 1 job ensures joblib does not copy arrays to sub-processes
    )
    
    dev_chunks.append(dev_chunk)
    
    del chunk_adata
    gc.collect()

# Concatenate chunks safely into a single AnnData object
print("Concatenating deviation chunks...")
final_chromvar_adata = sc.concat(dev_chunks, axis=0)

# Clean up chunk memory before adding to MuData
del dev_chunks
gc.collect()

# Add chromVAR matrix to mudata
mdata.mod['chromvar'] = final_chromvar_adata
# Add parameters
for parameter in [snakemake.params.sample_key, 'celltype']:
    barcode2param = mdata.obs[parameter].to_dict()
    mdata.mod['chromvar'].obs[parameter] = [barcode2param[x] for x in mdata.mod['chromvar'].obs_names]

mdata.write(snakemake.output.merged_multiome, compression='gzip')
