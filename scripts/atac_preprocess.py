import snapatac2 as snap

# Import and create AnnData object from fragment file
adata = snap.pp.import_fragments(
        snakemake.input.fragment_file, 
        file=None, 
        chrom_sizes=snap.genome.hg38.chrom_sizes,
        sorted_by_barcode=False,
        min_num_fragments=snakemake.params.min_peak_counts,
        n_jobs=snakemake.threads)

# Get the transcription start sites 
snap.metrics.tsse(adata, snap.genome.hg38)
snap.pp.filter_cells(adata, min_tsse=snakemake.params.min_tsse)

consensus_bed = snakemake.params.consensus_bed
if consensus_bed == None:
    snap.pp.add_tile_matrix(adata, bin_size=500)
else:
    snap.pp.add_peak_matrix(
        adata,
        peak_file = consensus_bed,
        inplace = True
    )

adata = adata[
        (adata.obs['n_fragment'] > snakemake.params.min_peak_counts) & 
        (adata.obs['tsse'] > snakemake.params.min_tsse)].copy()

# As this is a read-write interface, the AnnData object can just be closed
adata.write(snakemake.output.atac_anndata)
