import numpy as np
import pandas as pd
import gc
import scanpy as sc
import pychromvar as pc
import muon as mu
from pyjaspar import jaspardb
import scipy.sparse as sp

# Load MuData object
mdata = mu.read(snakemake.input.merged_multiome)

# CRITICAL: Keep memory light using float32 instead of float64
mdata['atac'].X = mdata['atac'].X.astype(np.float32)

# Really annoying script for running scipy version 1.17.1 with PyPi version of pychromVAR
def patch_scipy_sparse_sum():
    # Target both CSR and CSC legacy matrix base sum classes
    original_sum = sp._compressed._cs_matrix.sum

    def custom_sum(self, axis=None, dtype=None, out=None, keepdims=False):
        # If keepdims is invoked, convert to numpy/array to compute safely
        if keepdims:
            # Re-route the math using a modern numpy or array approach
            res = original_sum(self, axis=axis, dtype=dtype, out=out)
            # Reconstruct the expected 'keepdims' dimension shape
            if axis is not None:
                # Cast result safely to numpy, adjust shape, and return
                res_np = np.array(res)
                shape = list(self.shape)
                if isinstance(axis, int):
                    shape[axis] = 1
                else:
                    for ax in axis:
                        shape[ax] = 1
                return res_np.reshape(shape)
        # Fallback to standard behavior if keepdims isn't True
        return original_sum(self, axis=axis, dtype=dtype, out=out)

    # Bind the safe wrapper back directly into SciPy's source module in RAM
    sp._compressed._cs_matrix.sum = custom_sum
patch_scipy_sparse_sum()

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

# Update MuData object and save
mdata.mod['chromvar'] = final_chromvar_adata
mdata.write(snakemake.output.merged_multiome, compression='gzip')
