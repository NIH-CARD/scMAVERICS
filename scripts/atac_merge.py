import anndata as ad
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import scipy
import os

fragment_files = snakemake.input.fragment_file
output_files = snakemake.output.output_files# Read in fragments
adatas = snap.pp.import_fragments(
        fragment_files, 
        file = output_files,
        chrom_sizes=snap.genome.hg38.chrom_sizes,
        sorted_by_barcode=False,
        min_num_fragments=1000
        )

snap.metrics.tsse(adatas, snap.genome.hg38)

snap.pp.filter_cells(adatas, min_tsse=2.5)

if snakemake.input.consensus_bed != None:
        snap.pp.make_peak_matrix(
        adatas,
        inplace = True,
        peak_file = snakemake.input.consensus_bed,
        file = None,
        summary_type = 'sum'
        )

anndataset = snap.AnnDataSet(
    adatas=[(str(f.filename).split('/02_')[-1].split('_anndata_filtered_atac.h5ad')[0], f) for f in adatas], 
    filename=work_dir+'/atlas/temp_filtered_anndata_atac.h5ad')

dataset = snap.AnnDataSet(adatas=adatas, filename=snakemake.output.temp_file)
anndataset.obs['n_fragment'] = anndataset.adatas.obs['n_fragment']
anndataset.obs['tsse'] = anndataset.adatas.obs['tsse']

adata = anndataset.to_adata()
# Consolidate and export straight to a single permanent file
adata.write_h5ad(snakemake.output.merged_atac_anndata, compression = 'gzip')

os.remove(work_dir+'/atlas/temp_filtered_anndata_atac.h5ad')