import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly
data_dir = '/data/CARD_AUX/users/wellerca/PFC-atlas-preprocessing/CELLRANGER/'
#data_dir = '/data/CARD_singlecell/Brain_atlas/PFC_Multiome/CELLRANGER/'
# Define the working directory, explictly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
metadata_table = work_dir+'/input/metadata.csv'
# Define where celltypes/cell marker gene 
gene_markers_file = work_dir+'/input/first_pass_genes.csv'

# Key for samples, required in aggregating while preserving sample info
sample_key = 'SampleID'

# Read in the list of batches and samples
batches = pd.read_csv(metadata_table)['batch'].tolist()
samples = pd.read_csv(metadata_table)[sample_key].tolist()

# Name of the disease parameter
disease_param = 'cohort'
# Define disease states
control = 'NABEC'
diseases = ['HBCC']

# Define the cell types to look for, from gene marker file
cell_types = pd.read_csv(gene_markers_file)['cell type']

# Define RNA thresholds
mito_percent_thresh = 15
ribo_percent_thresh = 10
doublet_thresh = 0.15
min_genes_per_cell = 250

# Define ATAC thresholds
min_peak_counts = 250
min_num_cell_by_counts = 10


"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'singlecell': 'envs/single_cell_gpu.sif',
    'scenicplus': 'envs/scenicplus.sif',
    'snapatac': 'envs/snapATAC2.sif',
    'decoupler': 'envs/decoupler.sif',
    'circe': 'envs/circe.sif',
    'atac_fragment': 'envs/atac_fragment.sif'
    }

rule all:
    input:
        doublet_figure = work_dir+'/figures/QC_doublet.svg'
"""celltype_atac = expand(
    work_dir+'/hmmr_test/celltypes/{cell_type}/atac_circe.h5ad',
    cell_type = cell_types
    ),
merged_atac_anndata = work_dir+'/atlas/04_modeled_hmmr_atac.h5ad',
cell_specific_regression = expand(
    work_dir+'/hmmr_test/celltypes/{cell_type}/peak_age_regression.csv',
    cell_type = cell_types
)"""

"""merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu',
output_DGE_data = expand(
    work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
    cell_type = cell_types,
    disease = diseases
    ),
output_DAR_data = expand(
    work_dir + '/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
    cell_type = cell_types,
    disease = diseases
    ),
merged_cistopic_object = work_dir + '/data/pycisTopic/merged_cistopic_object.pkl',
merged_cistopic_adata = work_dir + '/atlas/05_annotated_anndata_atac.h5ad',
 rna_anndata=expand(
            data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            ),
atac_anndata = expand(
    data_dir+'{sample}/03_{sample}_anndata_object_atac.h5ad',
    zip,
    sample=samples,
    batch=batches
    ),"""
        
# This needs to be forced to run once
"""rule cellbender:
    input:
        rna_anndata =data_dir+'{sample}/raw_feature_bc_matrix.h5',
        cwd = data_dir+'{sample}'
    output:
        rna_anndata =data_dir+'{sample}/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=2880, mem_mb=300000, gpu=1, gpu_model='v100x'
    shell:
        work_dir+'/scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'

rule rna_preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'{sample}/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        sample='{sample}',
        sample_key = sample_key
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/rna_preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata=expand(
            data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'"""

rule plot_qc_rna:
    input:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    output:
        mito_figure = work_dir+'/figures/QC_mito_pct.svg',
        ribo_figure = work_dir+'/figures/QC_ribo_pct.svg',
        gene_counts_figure = work_dir+'/figures/QC_gene_counts.svg',
        doublet_figure = work_dir+'/figures/QC_doublet.svg',
        genes_by_counts = work_dir+'figures/QC_genes_by_counts.svg'
    singularity:
        envs['singlecell']
    resources:
        runtime=2880, mem_mb=3000000, disk_mb=10000, slurm_partition='largemem' 
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh,
        sample_key=sample_key,
    script:
        work_dir+'/scripts/plot_qc_metrics.py'

rule filter_rna:
    input:        
        rna_anndata =data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'{sample}/02_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
    script: 
        work_dir+'/scripts/rna_filter.py'

rule merge_filtered_rna:
    input:
        rna_anndata=expand(
            data_dir+'{sample}/02_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule atac_preprocess:
    input:
        fragment_file = data_dir+'{sample}/atac_fragments.tsv.gz'
    output:
        atac_anndata = data_dir+'{sample}/01_{sample}_anndata_object_atac.h5ad'
    singularity:
        envs['snapatac']
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/atac_preprocess.py'

rule merge_unfiltered_atac:
    input:
        rna_anndata=expand(
            data_dir+'{sample}/01_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_atac_anndata = work_dir+'/atlas/01_merged_anndata_atac.h5ad'
    singularity:
        envs['snapatac']
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_atac.py'

rule plot_qc_atac:
    input:
        atac_anndata=data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['snapatac']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/atac_plot_qc.py'

rule filter_atac:
    input:
        rna_anndata = data_dir+'{sample}/02_{sample}_anndata_filtered_rna.h5ad',
        atac_anndata = data_dir+'{sample}/01_{sample}_anndata_object_atac.h5ad'
    output:
        atac_anndata = data_dir+'{sample}/03_{sample}_anndata_object_atac.h5ad',
        rna_anndata = data_dir+'{sample}/03_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['snapatac']
    resources:
        runtime=30, mem_mb=50000, slurm_partition='quick'
    script:
        work_dir+'/scripts/atac_filter.py'

rule merge_multiome_rna:
    input:
        rna_anndata=expand(
            data_dir+'{sample}/03_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=480, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule feature_selection:
    input:
        merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
    output:
        hvg_rna_anndata = work_dir+'/atlas/03_hvg_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=360, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/feature_selection.py'

rule rna_model:
    input:
        hvg_rna_anndata = work_dir+'/atlas/03_hvg_anndata_rna.h5ad'
    output:
        hvg_rna_anndata = work_dir+'/atlas/04_modeled_hvg_anndata_rna.h5ad',
        model_history = work_dir+'/data/model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'/data/models/rna/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=4, gpu_model='v100x'
    shell:
        'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model}'

rule UMAP:
    input:
        merged_rna_anndata = work_dir + '/atlas/03_filtered_anndata_rna.h5ad',
        hvg_rna_anndata = work_dir + '/atlas/04_modeled_hvg_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir + '/atlas/04_modeled_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=1440, mem_mb=1000000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/scVI_to_UMAP.py'

rule first_pass_annotate:
    input:
        merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
        gene_markers = work_dir+'/input/first_pass_genes.csv'
    output:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
        cell_annotate = work_dir+'/data/first_pass_genes.csv'
    singularity:
        envs['singlecell']
    resources:
        runtime=480, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/annotate.py'

rule cluster_based_QC:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'/atlas/05_polished_anndata_rna.h5ad',
        course_celltype = work_dir + '/figures/first_pass_RNA_UMAP_celltype.svg',
        course_counts = work_dir + '/figures/first_pass_RNA_num_genes_celltype.svg'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir + '/scripts/cluster_based_QC.py'

rule filtered_feature_selection:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_polished_anndata_rna.h5ad'
    output:
        hvg_rna_anndata = work_dir+'/atlas/05_hvg_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=360, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/feature_selection.py'

rule rna_polish_model:
    input:
        hvg_rna_anndata = work_dir+'/atlas/05_hvg_anndata_rna.h5ad'
    output:
        hvg_rna_anndata = work_dir+'/atlas/05_modeled_hvg_anndata_rna.h5ad',
        model_history = work_dir+'/data/model_elbo/rna_model_v2_history.csv'
    params:
        model = work_dir+'/data/models/rna_V2/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=4, gpu_model='v100x'
    shell:
        'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model}'

rule filtered_UMAP:
    input:
        merged_rna_anndata = work_dir + '/atlas/05_polished_anndata_rna.h5ad',
        hvg_rna_anndata = work_dir + '/atlas/05_modeled_hvg_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir + '/atlas/06_modeled_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=1440, mem_mb=1000000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/scVI_to_UMAP.py'

rule second_pass_annotate:
    input:
        merged_rna_anndata = work_dir+'/atlas/06_modeled_anndata_rna.h5ad',
        gene_markers = gene_markers_file
    output:
        merged_rna_anndata = work_dir+'/atlas/07_annotated_anndata_rna.h5ad',
        cell_annotate = work_dir+'/data/example_marker_genes.csv'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/annotate.py'

rule gene_linear_regression:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_annotated_anndata_rna.h5ad',
        covariates = work_dir+'/data/covariates.csv'
    output:
        rna_pseudobulk = work_dir+'/data/pseudobulked_rna.csv',
        cell_gene_regression = work_dir+'/data/gene_age_regression.csv'
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/linear_regression_genes.py'

rule cistopic_pseudobulk:
    input:
        merged_rna_anndata = work_dir+'/atlas/06_modeled_anndata_rna.h5ad',
        fragment_file=expand(
            data_dir+'{sample}/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        pseudo_fragment_files = expand(
            work_dir + '/data/celltypes/{cell_types}/fragments.bed',
            cell_types=cell_types)
    params:
        pseudobulk_param = 'cell_type',
        samples=samples,
        sample_param_name = sample_key,
        cell_types = cell_types
    singularity:
        envs['atac_fragment']
    threads:
        64
    resources:
        runtime=960, mem_mb=3000000, disk_mb=500000, slurm_partition='largemem'
    script:
        'scripts/cistopic_pseudobulk.py'

rule MACS2_peak_call:
    input:
        pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/fragments.bed',
        blacklist_bed = work_dir + '/data/hg38-blacklist.v2.bed.gz'
    output: 
        narrow_peak = work_dir + "/hmmr_test/celltypes/{cell_type}/{cell_type}_accessible_regions.gappedPeak"
    params:
        out_dir = work_dir + "/data/celltypes/{cell_type}"
    resources:
        mem_mb=200000, runtime=2880
    shell:
        "macs2 callpeak --treatment {input.pseudo_fragment_files} --name {wildcards.cell_type} --outdir {params.out_dir} --format BEDPE --gsize hs --qvalue 0.001 --nomodel --shift 73 --extsize 146 --keep-dup all"

"module load macs/3; macs3 hmmratac --input {input.pseudo_fragment_files} --name {wildcards.cell_type} --outdir {params.out_dir} --format BEDPE "
"TEST CUTOFF ANALYSIS hmmratac --input ../fragments.bed --name astro --outdir . --blacklist ../../../hg38-blacklist.v2.bed.gz --format BEDPE --cutoff-analysis-only"

rule consensus_peaks:
    input:
        narrow_peaks = expand(
            work_dir + "/data/celltypes/{cell_type}/{cell_type}_accessible_regions.gappedPeak",
            cell_type = cell_types
            )
    output:
        consensus_bed = work_dir + '/data/consensus_regions.bed'
    singularity:
        envs['scenicplus']
    resources:
        runtime=960, mem_mb=100000
    script:
        'scripts/MACS_consensus.py'
        
rule cistopic_create_objects:
    input:
        merged_rna_anndata = work_dir+'/atlas/06_modeled_anndata_rna.h5ad',
        fragment_file = data_dir+'{sample}/atac_fragments.tsv.gz',
        consensus_bed = work_dir + '/data/consensus_regions.bed'
    output:
        cistopic_adata = data_dir + '{sample}/04_{sample}_anndata_peaks_atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample='{sample}'
    resources:
        runtime=120, mem_mb=350000
    threads:
        16
    script:
        'scripts/cistopic_create_object.py'

rule cistopic_merge_objects:
    input:
        atac_anndata=expand(
            data_dir + '{sample}/04_{sample}_anndata_peaks_atac.h5ad', 
            zip,
            sample=samples,
            batch=batches
            )
    output:
        merged_atac_anndata = work_dir + '/atlas/03_merged_hmmr_atac.h5ad'
    singularity:
        envs['atac_fragment']
    params:
        sample_key = sample_key,
        samples = samples
    resources:
        runtime=2880, mem_mb=1000000, slurm_partition='largemem'
    script:
        'scripts/merge_cistopic_and_adata.py'

rule atac_peaks_model:
    input:
        merged_atac_anndata = work_dir+'/atlas/03_merged_cistopic_atac.h5ad'
    output:
        merged_atac_anndata = work_dir+'/atlas/04_modeled_cistopic_atac.h5ad',
        atac_model_history = work_dir+'/data/model_elbo/atac_cistopic_model_history.csv'
    params:
        atac_model = work_dir+'/data/models/atac/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/atac_model.sh {input.merged_atac_anndata} {params.sample_key} {output.atac_model_history} {output.merged_atac_anndata} {params.atac_model}'

rule multiome_output:
    input:
        merged_atac_anndata = work_dir + '/atlas/04_modeled_cistopic_atac.h5ad',
        merged_rna_anndata = work_dir+'/atlas/06_modeled_anndata_rna.h5ad'
    output:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    singularity:
        envs['singlecell']
    resources:
        runtime=120, mem_mb=300000, slurm_partition='quick' 
    script:
        'scripts/merge_muon.py'

rule create_bigwig:
    input:
        pseudo_fragment_file = work_dir + '/data/celltypes/{cell_type}/fragments.bed'
    output:
        celltype_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_bigwig.bw',
        celltype_normalized_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_normalized_bigwig.bw'
    resources:
        mem_mb=1500000, runtime=960, slurm_partition='largemem'
    singularity:
        envs['atac_fragment']
    script:
        'scripts/atac_bigwig.py'

rule celltype_bed:
    input:
        xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_accessible_regions.gappedPeak",
    output:
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed'
    singularity:
        envs['atac_fragment']
    script:
        'scripts/MACS_to_bed.py'

rule annotate_bed:
    input:
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed'
    output:
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed'
    resources:
        runtime=30, mem_mb=50000,  slurm_partition='quick,norm'
    shell:
        'module load homer;annotatePeaks.pl {input.cell_bedfile} hg38 > {output.cell_annotated_bedfile}'

rule export_atac_cell:
    input:
        merged_rna_anndata = work_dir+'/atlas/06_modeled_anndata_rna.h5ad',
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed',
        fragment_files=expand(
            data_dir+'{sample}/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        samples=samples,
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    threads:
        8
    resources:
        runtime=2880, mem_mb=800000, slurm_partition='largemem'
    script:
        'scripts/atac_by_celltype.py'


rule peak_linear_regression:
    input:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad',
        covariates = '/data/CARD_singlecell/PFC_atlas/data/covariates.csv'
    output:
        cell_specific_pseudo = work_dir+'/data/celltypes/{cell_type}/pseudobulk_atac.csv',
        cell_specific_regression = work_dir+'/data/celltypes/{cell_type}/peak_age_regression.csv'
    params:
        cell_type = lambda wildcards, output: output[0].split("_")[-2],
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/linear_regression_peaks.py'

rule wnn:
    input:
        merged_atac_anndata = work_dir + '/atlas/04_modeled_cistopic_atac.h5ad',
        merged_rna_anndata = work_dir+'/atlas/06_modeled_anndata_rna.h5ad'
    output:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    singularity:
        envs['snapatac']
    resources:
        runtime=960, mem_mb=1500000, slurm_partition='largemem' 
    script:
        'scripts/wnn.py'

rule atac_coaccessibilty:
    input:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
    output:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac_circe.h5ad'
    params:
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    singularity:
        envs['circe']
    threads:
        8
    resources:
        runtime=5760, mem_mb=1000000, slurm_partition='largemem'
    script:
        'scripts/circe_by_celltype.py'