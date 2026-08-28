
import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

configfile: "config.yaml"
"""File locations"""
data_dir = '/data/CARD_singlecell/Brain_atlas/PCA_Multiome/' # Define the data directory, explicitly
work_dir = '/data/CARD_singlecell/PCA_multiome/' # Define the working directory, explictly as the directory of this pipeline
metadata_table = work_dir+'input/sequenced_samples.csv' # Define where the metadata data exists for each sample to be processed
gene_markers_file = work_dir+'input/first_pass_genes.csv' # Define where celltypes/cell marker gene 

"""Metadata parameters"""
seq_batch_key = 'batch' # Key for sequencing batch, used for directory search`
sample_key = 'CARD_ID' # Key for samples, required in aggregating while preserving sample info


samples = pd.read_csv(metadata_table)[sample_key].tolist()
disease_param = 'Pathology' # Name of the disease parameter
control = 'control' # Define disease states
diseases = ['PD', 'DLB'] # Disease states to compare, keep as list of strings, unnecessary 
disease_comparisons = ['control vs. PD', 'control vs. DLB', 'PD vs. DLB']
cell_types = pd.read_csv(gene_markers_file)['cell type'] # Define the cell types to look for, from gene marker file
design_covariates = ['Age','Sex'] # Design factors/covariates for DGEs and DARs
reference_genome = '/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa' 
genome_length = '/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/star/chrNameLength.txt'
cell_cycle_gene_file = work_dir + 'input/lab_cell_cycle_genes.txt'

"""Quality control thresholds"""
mito_percent_thresh = 15 # Maximum percent of genes in a cell that can be mitochondrial
ribo_percent_thresh = 10 # Maximum percent of genes in a cell that can be ribosomal
doublet_thresh = 0.15 # Maximum doublet score for a cell, computed by scrublet
min_genes_per_cell = 250 # Minimum number of unique genes in a cell
min_peak_counts = 1000 # Minimum number of fragments per cell
min_tsse = 2.5 # Minimum transcription start site enrichment

""" Samples processed so far, remove once all samples have been sequenced"""
working_samples = pd.read_csv(work_dir + 'input/sequenced_samples.csv')['CARD_ID'].to_list()
working_batches = pd.read_csv(work_dir + 'input/sequenced_samples.csv')['batch'].to_list()

batches = working_batches # Read in the list of batches and samples
samples = working_samples
subtypes = []

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'snapatac2': 'envs/snapatac2.sif',
    'singlecell': 'envs/single_cell_gpu.sif',
    'tobias': 'envs/tobias.sif',
    'dreampy': 'envs/dreampy.sif',
    'multiome': 'envs/multiome.sif',
    'scenicplus': 'envs/scenicplus.sif',
    }

rule all:
    input:
        merged_rna_anndata = work_dir+'atlas/05_QC_filtered_anndata_rna.h5ad'

# This needs to be forced to run once
"""rule cellbender:
    input:
        rna_anndata =data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/raw_feature_bc_matrix.h5',
        cwd = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/'
    output:
        rna_anndata = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=1440, mem_mb=64000, gpu=1, gpu_model='v100x'
    shell:
        work_dir+'/scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'
"""

"""========================================================================="""
"""                                RNA portion                              """
"""========================================================================="""

"""
rule cellbender:
    input:
        rna_anndata =data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/raw_feature_bc_matrix.h5',
        cwd = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/'
    output:
        rna_anndata = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=1440, mem_mb=200000, gpu=1, gpu_model='v100x'
    shell:
        work_dir+'scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'
"""
rule rna_preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['multiome']
    params:
        sample='{sample}',
        sample_key = sample_key,
        cell_cycle_gene_file = cell_cycle_gene_file
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'scripts/rna_preprocess.py'

rule rna_filter:
    input:        
        rna_anndata = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['multiome']
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, slurm_partition='quick' 
    script: 
        work_dir+'scripts/rna_filter.py'

rule rna_merge:
    input:
        rna_anndata=expand(
            data_dir+'batch{batch}/cellranger/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'atlas/02_filtered_anndata_rna.h5ad'
    singularity:
        envs['multiome']
    params:
        samples=samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/rna_merge.py'

rule rna_feature_selection:
    input:
        merged_rna_anndata = work_dir+'atlas/02_filtered_anndata_rna.h5ad'
    output:
        hvg_rna_anndata = work_dir+'atlas/03_hvg_anndata_rna.h5ad'
    singularity:
        envs['multiome']
    params:
        num_hvgenes = config['high_var_gene_num']
    resources:
        runtime=360, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'scripts/rna_feature_selection.py'

rule rna_model:
    input:
        hvg_rna_anndata = work_dir+'atlas/03_hvg_anndata_rna.h5ad'
    output:
        hvg_rna_anndata = work_dir+'atlas/04_modeled_hvg_anndata_rna.h5ad',
        model_history = work_dir+'data/model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'data/models/rna/',
        sample_key = sample_key,
        random_number_seed = config['random_number_seed'],
        num_layers = config['num_layers'],
        num_latent = config['num_latent'],
        max_epoch = config['max_epoch'],
        machine_type = config['machine_type']
    resources:
        runtime=2880, mem_mb=200000, gpu=2, gpu_model="a100"
    shell:
        'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model} {params.random_number_seed} {params.num_layers} {params.num_latent} {params.max_epoch} {params.machine_type}'

rule rna_latent_transfer:
    input:
        merged_rna_anndata = work_dir + 'atlas/02_filtered_anndata_rna.h5ad',
        hvg_rna_anndata = work_dir + 'atlas/04_modeled_hvg_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir + 'atlas/04_modeled_anndata_rna.h5ad'
    singularity:
        envs['multiome']
    resources:
        runtime=1440, mem_mb=1000000, slurm_partition='largemem'
    script:
        work_dir+'scripts/rna_latent_transfer.py'

rule rna_annotate:
    input:
        merged_rna_anndata = work_dir+'atlas/04_modeled_anndata_rna.h5ad',
        gene_markers = work_dir+'input/first_pass_genes.csv'
    output:
        merged_rna_anndata = work_dir+'atlas/05_annotated_anndata_rna.h5ad',
    params:
        seq_batch_key = seq_batch_key
    singularity:
        envs['multiome']
    resources:
        runtime=480, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'scripts/rna_annotate.py'

rule rna_cluster_based_QC:
    input:
        merged_rna_anndata = work_dir+'atlas/05_annotated_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'atlas/05_QC_filtered_anndata_rna.h5ad',
        course_celltype = work_dir + 'figures/first_pass_RNA_UMAP_celltype.svg',
        course_counts = work_dir + 'figures/first_pass_RNA_num_genes_celltype.svg'
    singularity:
        envs['multiome']
    resources:
        runtime=240, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir + '/scripts/rna_cluster_based_QC.py'