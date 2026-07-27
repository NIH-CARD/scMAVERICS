import snapatac2 as snap

# Import and create AnnData object from fragment file
adata = snap.pp.import_data(
        snakemake.input.fragment_file, 
        file=None, 
        chrom_sizes=snap.genome.hg38.chrom_sizes,
        sorted_by_barcode=False,
        min_num_fragments=1,
        n_jobs=snakemake.threads)

# Get the transcription start sites 
snap.metrics.tsse(adata, snap.genome.hg38)

snap.pp.add_tile_matrix(adata, bin_size=500)
snap.pp.select_features(adata, n_features=50000)
snap.pp.scrublet(adata, n_comps=5, sim_doublet_ratio=0.5)
# As this is a read-write interface, the AnnData object can just be closed
adata.write(snakemake.output.atac_anndata, compression='gzip')
