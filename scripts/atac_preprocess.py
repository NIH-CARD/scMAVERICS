import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as snap
import pandas as pd

# Import and create AnnData object from fragment file
adata = snap.pp.import_data(
        snakemake.input.fragment_file, 
        file=snakemake.output.atac_anndata, 
        chrom_sizes=snap.genome.hg38.chrom_sizes,
        sorted_by_barcode=False,
        n_jobs=2)

# Get the fragment distribution (for later QC)
#_ = snap.pl.frag_size_distr(adata, show=False)
# Get the transcription start sites 
snap.metrics.tsse(adata, snap.genome.hg38)

# Add tiled bins across chromosomes, bin size of 500
snap.pp.add_tile_matrix(adata, bin_size=500)
# For this sample, determine the count and highly variable bins
snap.pp.select_features(adata, n_features=250000)
# As this is a read-write interface, the AnnData object can just be closed
adata.close()