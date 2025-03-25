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
metadata_table = work_dir+'/input/example_metadata.csv'
# Define where celltypes/cell marker gene 
gene_markers_file = work_dir+'/input/example_marker_genes.csv'

# Key for samples, required in aggregating while preserving sample info
sample_key = 'Sample'

# Read in the list of batches and samples
batches = pd.read_csv(metadata_table)['Use_batch'].tolist()
samples = pd.read_csv(metadata_table)[sample_key].tolist()

# Name of the disease parameter
disease_param = 'Primary Diagnosis'
# Define disease states
control = 'control'
diseases = ['PD', 'DLB']

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
    'singlecell': 'envs/single_cell_cpu.sif', 
    'single_cell_gpu': 'envs/single_cell_gpu.sif',
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
        merged_cistopic_adata = work_dir + '/atlas/05_annotated_anndata_atac.h5ad',
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
        
rule cellbender:
    input:
        rna_anndata =data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/raw_feature_bc_matrix.h5',
        cwd = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=2880, mem_mb=300000, gpu=1, gpu_model='v100x'
    shell:
        work_dir+'/scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'

rule rna_preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        sample='{sample}'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/rna_preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad', 
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
    singularity:
        envs['singlecell']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/plot_qc_metrics.py'

rule filter_rna:
    input:        
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad'
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
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad', 
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
        fragment_file=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz'
    output:
        atac_anndata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/atac_preprocess.py'

rule merge_unfiltered_atac:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_atac_anndata = work_dir+'/atlas/01_merged_anndata_atac.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_atac.py'

rule plot_qc_atac:
    input:
        atac_anndata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/atac_plot_qc.py'

rule filter_atac:
    input:
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad',
        atac_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad'
    output:
        atac_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad',
        rna_anndata = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=30, mem_mb=50000, slurm_partition='quick'
    script:
        work_dir+'/scripts/atac_filter.py'

rule merge_multiome_rna:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad', 
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

rule rna_model:
    input:
        merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
        model_history = work_dir+'/model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'/data/models/rna/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/rna_model.sh {input.merged_rna_anndata} {params.sample_key} {output.model_history} {output.merged_rna_anndata} {params.model}'

rule annotate:
    input:
        merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
        gene_markers = gene_markers_file
    output:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
        cell_annotate = work_dir+'/data/rna_cell_annot.csv'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'

"""LEGACY BIN METHOD"""

"""rule atac_bins_model:
    input:
        cell_annotate = work_dir+'/data/rna_cell_annot.csv',
        metadata_table=metadata_table,
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        umap_data = work_dir+'/data/atac_umap.csv',
        var_data = work_dir+'/data/atac_var_selected.csv',
        merged_atac_anndata = work_dir+'/atlas/03_filtered_anndata_atac.h5ad'
    params:
        samples = samples,
        sample_key = sample_key
    singularity:
        envs['singlecell']
    resources:
        runtime=2880, mem_mb=3000000, slurm_partition='largemem'
    threads:
        64
    script:
        work_dir+'/scripts/merge_atac.py' 

rule atac_bins_annotate:
    input:
        atac_anndata = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad', 
            zip,
            batch=batches,
            sample=samples
            ),
        umap_csv = work_dir+'/data/atac_umap.csv',
        var_csv = work_dir+'/data/atac_var_selected.csv',
        annot_csv = work_dir+'/data/rna_cell_annot.csv',
        metadata_table = metadata_table
    output:
        atac_anndata = work_dir+'/atlas/04_filtered_anndata_atac.h5ad',
        merged_atac_anndata = work_dir+'/atlas/05_annotated_anndata_atac.h5ad'
    params:
        samples=samples
    singularity:
        envs['singlecell']
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000
    script:
        'scripts/atac_bins_annotate.py'"""

rule DGE:
    input:
        rna_anndata = work_dir + '/atlas/05_annotated_anndata_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
        output_figure = work_dir + 'figures/{cell_type}/rna_{cell_type}_{disease}_DAR.png',
        celltype_pseudobulk = work_dir+'/data/celltypes/{cell_type}/rna_{cell_type}_{disease}_pseudobulk.h5ad'
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
        'scripts/rna_DGE.py'

rule cistopic_pseudobulk:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
        fragment_file=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        bigwig_paths = work_dir + '/data/pycisTopic/pseudobulk_bigwig_files/bw_paths.tsv',
        bed_paths = work_dir + '/data/pycisTopic/pseudobulk_bed_files/bed_paths.tsv'
    params:
        bigwig_file_locs = work_dir + '/data/pycisTopic/pseudobulk_cell_bigwig_files/',
        bed_file_locs = work_dir + '/data/pycisTopic/pseudobulk_cell_bed_files/',
        pseudobulk_param = 'cell_type',
        samples=samples,
        sample_param_name = sample_key
    singularity:
        envs['scenicplus']
    threads:
        64
    resources:
        runtime=240, mem_mb=3000000, disk_mb=500000, slurm_partition='largemem'
    script:
        'scripts/cistopic_pseudobulk.py'

rule cistopic_call_peaks:
    input:
        bigwig_paths = work_dir + '/data/pycisTopic/pseudobulk_bigwig_files/bw_paths.tsv',
        bed_paths = work_dir + '/data/pycisTopic/pseudobulk_bed_files/bed_paths.tsv'
    output:
        consensus_bed = work_dir + '/data/pycisTopic/consensus_regions.bed',
        peak_dict = work_dir + '/data/pycisTopic/MACS/narrow_peaks_dict.pkl'
    params:
        MACS_dir = work_dir + '/data/pycisTopic/MACS'
    singularity:
        envs['scenicplus']
    resources:
        runtime=240, mem_mb=100000, disk_mb=500000
    script:
        'scripts/cistopic_call_peaks.py'
    
rule cistopic_create_objects:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
        fragment_file = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
        consensus_bed = work_dir + '/data/pycisTopic/consensus_regions.bed'
    output:
        cistopic_object = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
        cistopic_adata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample='{sample}'
    resources:
        runtime=240, mem_mb=200000
    threads:
        16
    script:
        'scripts/cistopic_create_object.py'

rule cistopic_merge_objects:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
        cistopic_objects = expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
            zip,
            sample=samples,
            batch=batches
            ),
        rna_anndata=expand(
            data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad', 
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
        runtime=1440, mem_mb=2000000, slurm_partition='largemem'
    script:
        'scripts/merge_cistopic_and_adata.py'

rule atac_peaks_model:
    input:
        merged_atac_anndata = work_dir+'/atlas/03_merged_cistopic_atac.h5ad'
    output:
        merged_atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad',
        atac_model_history = work_dir+'/model_elbo/atac_model_history.csv'
    params:
        atac_model = work_dir+'/data/models/atac/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/atac_model.sh {input.merged_rna_anndata} {params.sample_key} {output.atac_model_history} {output.merged_rna_anndata} {params.atac_model}'

rule atac_peaks_annotate:
    input:
        merged_atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad',
        annot_csv = work_dir+'/data/rna_cell_annot.csv'
    output:
        merged_cistopic_adata = work_dir + '/atlas/05_annotated_anndata_atac.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=300000
    script:
        'scripts/atac_annotate.py'

rule DAR:
    input:
        atac_anndata = work_dir+'data/celltypes/{cell_type}/atac.h5ad'
    output:
        output_DAR_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
        output_figure = work_dir+'/figures/{cell_type}/atac_{cell_type}_{disease}_DAR.png'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['singlecell']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/atac_DAR.py'
   
rule multiome_output:
    input:
        merged_cistopic_adata = work_dir + '/atlas/05_annotated_anndata_atac.h5ad',
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad'
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
        envs['singlecell']
    threads:
        8
    resources:
        runtime=120, mem_mb=300000
    script:
        'scripts/export_celltype.py'

