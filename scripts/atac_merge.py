import anndata as ad
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import scipy

# Read in anndata objects using backed mode
adata_files = snakemake.output.output_files
fragments = snakemake.input.fragment_file
sample_names = snakemake.params.samples
consensus_bed = snakemake.input.consensus_bed

adatas = snap.pp.import_fragments(
        fragments, 
        file=adata_files, 
        chrom_sizes=snap.genome.hg38.chrom_sizes,
        sorted_by_barcode=False,
        min_num_fragments=1000,
        n_jobs = snakemake.threads
        )

# Transcription start site filter
snap.metrics.tsse(adatas, snap.genome.hg38, n_jobs=snakemake.threads)
# Minimum number of fragments filter
snap.pp.filter_cells(adatas, min_tsse=2.5, n_jobs = snakemake.threads)

adataset = snap.AnnDataSet(
    adatas=[(sample_names[i], adatas[i]) for i in (rangelen(adatas))],
    filename=snakemake.output.temp_merged_anndata
)
print('Dataset merged')
adataset.obs_names = np.array(adataset.obs_names) + '_' + np.array(adataset.obs['sample'])
adataset.obs['n_fragment'] = adataset.adatas.obs['n_fragment']
adataset.obs['tsse'] = adataset.adatas.obs['tsse']

# Populate merged object with calculated metadata
if consensus_bed != None:
    atac = snap.pp.make_peak_matrix(
        adataset,
        peak_file = consensus_bed
        )

# Close merged output file
atac.write_h5ad(snakemake.output.merged_atac_anndata, compression = 'gzip')
print('Merged anndata object closed')