# Import general packages
import numpy as np
import pandas as pd

# Import specialized packages
import scanpy as sc
import pychromvar as pc
import muon as mu
from pyjaspar import jaspardb

# Load MuData object
mdata = mu.read(snakemake.input.merged_multiome)

# Get reference genome
pc.add_peak_seq(
    mdata,
    genome_file = snakemake.input.reference_genome,
    delimiter = ':|-'
)

# Correct for GC bias
pc.add_gc_bias(mdata)
pc.get_bg_peaks(
    mdata,
    n_jobs = snakemake.threads
)

# Fetch motifs (modify if non-human)
jdb_obj = jaspardb(release='JASPAR2024')
motifs = jdb_obj.fetch_motifs(
    collection = 'CORE',
    tax_group = ['vertebrates'])

# Match motifs in each sample
pc.match_motif(mdata, motifs=motifs)

# Return expectation values 
dev = pc.compute_deviations(mdata, n_jobs = -1)

# Save the motif value AnnData inside the muon object
mdata.mod['chromvar'] = dev
mdata.mod['chromvar'].raw = dev

# Write the muon object
mdata.write(snakemake.output.merged_multiome, compression='gzip')