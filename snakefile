import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly
data_dir = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/'
# Define the working directory, explictly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
metadata_table = work_dir+'/input/SN_PD_DLB_samples.csv'
# Define where celltypes/cell marker gene 
gene_markers_file = work_dir+'/input/example_marker_genes.csv'

# Key for samples, required in aggregating while preserving sample info
sample_key = 'Sample'

# Read in the list of batches and samples
batches = pd.read_csv(metadata_table)['Use_batch'].tolist()
samples = pd.read_csv(metadata_table)[sample_key].tolist()
disease_param = 'Primary Diagnosis' # Name of the disease parameter
control = 'control' # Define disease states
diseases = ['PD', 'DLB'] # Disease states to compare, keep as list of strings, unnecessary 
cell_types = pd.read_csv(gene_markers_file)['cell type'] # Define the cell types to look for, from gene marker file
design_covariates = [seq_batch_key, 'Age', 'Sex'] # Design factors/covariates for DGEs and DARs

"""Quality control thresholds"""
mito_percent_thresh = 15 # Maximum percent of genes in a cell that can be mitochondrial
ribo_percent_thresh = 10 # Maximum percent of genes in a cell that can be ribosomal
doublet_thresh = 0.15 # Maximum doublet score for a cell, computed by scrublet
min_genes_per_cell = 250 # Minimum number of unique genes in a cell
min_peak_counts = 500 # Minimum number of fragments per cell


"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'snapatac2': 'envs/snapatac2.sif',
    'singlecell': 'envs/single_cell_gpu.sif',
    'scenicplus': 'envs/scenicplus.sif',
    'decoupler': 'envs/decoupler.sif'
    }

rule all:
    input:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu',
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
        merged_cistopic_adata = work_dir + '/atlas/05_annotated_cistopic_atac.h5ad',
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            ),
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad',
            zip,
            sample=samples,
            batch=batches
            ),
        
# This needs to be forced to run once
rule cellbender:
    script:
        work_dir+'/scripts/cellbender_array.sh'

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
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule plot_qc_rna:
    input:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    output:
        mito_figure = work_dir+'/figures/QC_mito_pct.png',
        ribo_figure = work_dir+'/figures/QC_ribo_pct.png',
        gene_counts_figure = work_dir+'/figures/QC_gene_counts.png',
        doublet_figure = work_dir+'/figures/QC_doublet.png',
        genes_by_counts = work_dir+'figures/QC_genes_by_counts.png'
    singularity:
        envs['singlecell']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
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
        rna_anndata = data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad'
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
        fragment_file=data_dir+'{sample}/atac_fragments.tsv.gz'
    output:
        atac_anndata=data_dir+'{sample}/01_{sample}_anndata_object_atac.h5ad'
    singularity:
        envs['snapatac2']
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
        envs['snapatac2']
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_atac.py'

rule plot_qc_atac:
    input:
        atac_anndata=data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['snapatac2']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    params:
        sample_key=sample_key
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
        envs['snapatac2']
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
        runtime=120, mem_mb=300000, disk_mb=10000#, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule feature_selection:
    input:
        merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
        model_history = work_dir+'/model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'/data/models/rna_v2/',
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
    params:
        seq_batch_key = seq_batch_key
    singularity:
        envs['singlecell']
    resources:
        runtime=480, mem_mb=1500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'

rule DGE:
    input:
        rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
        output_figure = work_dir + '/figures/{cell_type}/rna_{cell_type}_{disease}_DGE.svg',
        celltype_pseudobulk = work_dir+'/data/celltypes/{cell_type}/rna_{cell_type}_{disease}_pseudobulk.csv'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3],
        sample_key=sample_key,
        design_factors = design_covariates,
        separating_cluster = 'cell_type'
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/rna_DGE.py'


rule cistopic_pseudobulk:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        fragment_file=expand(
            data_dir+'{sample}/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        bigwig_paths = work_dir + '/data/pycisTopic/pseudobulk_bigwig_files/bw_paths.tsv',
        bed_paths = work_dir + '/data/pycisTopic/pseudobulk_bed_files/bed_paths.tsv'
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
        runtime=240, mem_mb=3000000, disk_mb=500000, slurm_partition='largemem'
    script:
        'scripts/cistopic_pseudobulk.py'

rule cistopic_call_peaks:
    input:
        pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/{cell_type}_fragments.bed'
    output: 
        xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.xls",
        narrow_peak = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.narrowPeak"
    params:
        out_dir = work_dir + "/data/celltypes/{cell_type}"
    resources:
        mem_mb=200000, runtime=960
    singularity:
        envs['scenicplus']
    shell:
        "macs2 callpeak --treatment {input.pseudo_fragment_files} --name {wildcards.cell_type} --outdir {params.out_dir} --format BEDPE --gsize hs --qvalue 0.001 --nomodel --shift 73 --extsize 146 --keep-dup all"

rule consensus_peaks:
    input:
        narrow_peaks = expand(
            work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.narrowPeak",
            cell_type = cell_types
            )
    output:
        consensus_bed = work_dir + '/data/consensus_regions.bed'
    singularity:
        envs['scenicplus']
    threads:
        32
    resources:
        runtime=240, mem_mb=100000, disk_mb=500000
    script:
        'scripts/MACS_consensus.py'
    
rule cistopic_create_objects:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        fragment_file = data_dir+'{sample}/atac_fragments.tsv.gz',
        consensus_bed = work_dir + '/data/consensus_regions.bed'
    output:
        cistopic_object = data_dir+'{sample}/04_{sample}_cistopic_obj.pkl',
        cistopic_adata = data_dir+'{sample}/04_{sample}_anndata_peaks_atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample='{sample}',
        seq_batch_key = seq_batch_key,
        sample_key = sample_key,
        disease_param = disease_param
    resources:
        runtime=120, mem_mb=250000, slurm_partition='quick'
    threads:
        16
    script:
        'scripts/cistopic_create_object.py'

rule cistopic_merge_objects:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        cistopic_objects = expand(
            data_dir+'{sample}/04_{sample}_cistopic_obj.pkl',
            zip,
            sample=samples,
            batch=batches
            ),
        rna_anndata=expand(
            data_dir+'{sample}/04_{sample}_anndata_peaks_atac.h5ad', 
            zip,
            sample=samples,
            batch=batches
            )
    output:
        merged_cistopic_object = work_dir + '/data/pycisTopic/merged_cistopic_object.pkl',
        merged_cistopic_adata = work_dir + '/atlas/03_merged_cistopic_atac.h5ad'
    singularity:
        envs['scenicplus']
    resources:
        runtime=960, mem_mb=300000
    script:
        'scripts/merge_cistopic_and_adata.py'

rule atac_peaks_model:
    input:
        merged_atac_anndata = work_dir+'/atlas/03_merged_cistopic_atac.h5ad'
    output:
        merged_atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad',
        atac_model_history = work_dir+'/data/model_elbo/atac_model_history.csv'
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
        merged_atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad',
        annot_csv = work_dir+'/data/rna_cell_annot.csv'
    output:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=300000
    script:
        'scripts/export_celltype.py'


rule DAR:
    input:
        atac_anndata = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
    output:
        output_DAR_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
        output_figure = work_dir+'/figures/{cell_type}/atac_{cell_type}_{disease}_DAR.svg',
        cell_specific_pseudo = work_dir+'/data/celltypes/{cell_type}/atac_{disease}_pseudobulk.csv'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/atac_DARs.py'
   
rule atac_coaccessibilty:
    input:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
    output:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    singularity:
        envs['singlecell']
    script:
        'scripts/merge_muon.py'

rule export_celltypes:
    input:
        merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
    output:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/atac.h5ad',
        celltype_rna = work_dir+'data/celltypes/{cell_type}/rna.h5ad'
    params:
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    singularity:
        envs['circe']
    threads:
        8
    resources:
        runtime=2880, mem_mb=1500000, slurm_partition='largemem'
    script:
        'scripts/circe_by_celltype.py'
