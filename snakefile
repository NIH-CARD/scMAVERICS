import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""


"""File locations"""
data_dir = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/' # Define the data directory, explicitly
work_dir = '/data/CARD_singlecell/SN_atlas' # Define the working directory, explictly as the directory of this pipeline
metadata_table = work_dir+'/input/SN_PD_DLB_samples.csv' # Define where the metadata data exists for each sample to be processed
gene_markers_file = work_dir+'/input/example_marker_genes.csv' # Define where celltypes/cell marker gene 

"""Metadata parameters"""
seq_batch_key = 'Use_batch' # Key for sequencing batch, used for directory search
sample_key = 'Sample_ID' # Key for samples, required in aggregating while preserving sample info
batches = pd.read_csv(metadata_table)[seq_batch_key].tolist() # Read in the list of batches and samples

samples = pd.read_csv(metadata_table)[sample_key].tolist()
disease_param = 'Primary Diagnosis' # Name of the disease parameter
control = 'control' # Define disease states
diseases = ['PD', 'DLB'] # Disease states to compare, keep as list of strings, unnecessary 
cell_types = pd.read_csv(gene_markers_file)['cell type'] # Define the cell types to look for, from gene marker file
design_covariates = ['Age','Sex'] # Design factors/covariates for DGEs and DARs
reference_genome = '/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa' 
genome_length = '/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/star/chrNameLength.txt'

"""Quality control thresholds"""
mito_percent_thresh = 15 # Maximum percent of genes in a cell that can be mitochondrial
ribo_percent_thresh = 10 # Maximum percent of genes in a cell that can be ribosomal
doublet_thresh = 0.15 # Maximum doublet score for a cell, computed by scrublet
min_genes_per_cell = 250 # Minimum number of unique genes in a cell
min_peak_counts = 500 # Minimum number of fragments per cell

"""Subcluster values, extracted manually after"""
leiden_clusters =['0', '6', '1', '5', '24', '26', '20', '13', '28', '32', '27', '34', '15', '11', '30', '7', '21', '36', '4', '19', '25', '37', '8', '17', '18', '9', '14', '38', '35', '3', '33', '31', '16', '39', '41', '42', '2', '23', '43', '44', '45','40']

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'snapatac2': 'envs/snapatac2.sif',
    'singlecell': 'envs/single_cell_gpu.sif',
    'scenicplus': 'envs/scenicplus.sif',
    'decoupler': 'envs/decoupler.sif',
    'circe': 'envs/circe.sif',
    'atac_fragment': 'envs/atac_fragment.sif',
    'great_gsea': 'envs/great_gsea.sif'
    }

rule all:
    input:
        footprinted_bigwig = expand(
            work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_footprints.bw',
            cell_type = cell_types,
            disease = ['PD', 'DLB', 'control']
        ),
        control_footprint_bigwig = expand(
            work_dir+'/data/celltypes/{cell_type}/{cell_type}_control_ATACorrect/{cell_type}_control_comparison_footprints.bw',
            cell_type = cell_types
        )
        
"""
output_DGE_data = expand(
    work_dir + '/data/significant_genes/rna/leiden/rna_{cell_type}_PD_vs_{disease}_DGE.csv',
    cell_type = leiden_clusters,
    disease = ['DLB']
    ),
output_leiden_DAR_data = expand(
    work_dir+'/data/significant_genes/atac/leiden/atac_{cell_type}_{disease}_DAR.csv',
    cell_type = leiden_clusters,
    disease = diseases
),"""
        

# This needs to be forced to run once
rule cellbender:
    input:
        rna_anndata =data_dir+'{sample}/raw_feature_bc_matrix.h5',
        cwd = data_dir+'{sample}/'
    output:
        rna_anndata = data_dir+'{sample}/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=1440, mem_mb=300000, gpu=2, gpu_model='v100x'
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
        work_dir+'/scripts/annotate.py'

rule cluster_based_QC:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'/atlas/05_QC_filtered_anndata_rna.h5ad',
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
        merged_rna_anndata = work_dir+'/atlas/05_QC_filtered_anndata_rna.h5ad'
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
        model = work_dir+'/data/models/rna_polish/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model}'

rule filtered_UMAP:
    input:
        merged_rna_anndata = work_dir + '/atlas/05_QC_filtered_anndata_rna.h5ad',
        hvg_rna_anndata = work_dir + '/atlas/05_modeled_hvg_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir + '/atlas/06_polished_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=1440, mem_mb=1000000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/scVI_to_UMAP.py'

"""rule second_pass_annotate:
    input:
        merged_rna_anndata = work_dir+'/atlas/06_polished_anndata_rna.h5ad',
        gene_markers = gene_markers_file
    output:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        cell_annotate = work_dir+'/data/rna_cell_annot.csv'
    params:
        seq_batch_key = seq_batch_key
    singularity:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/annotate.py'"""

rule gene_linear_regression:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        covariates = work_dir+'/data/covariates.csv'
    output:
        rna_pseudobulk = work_dir+'/data/pseudobulked_rna.csv',
        cell_gene_regression = work_dir+'/data/gene_age_regression.csv'
    params:
        sample_key=sample_key,
        disease_param = disease_param,
        design_factors = design_covariates,
        cell_types = cell_types
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/linear_regression_genes.py'

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
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        pseudo_fragment_files = expand(
            work_dir + '/data/celltypes/{cell_type}/{cell_type}_fragments.bed',
            cell_type = cell_types)
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
        'scripts/fragment_pseudobulk.py'

rule cistopic_call_peaks:
    input:
        pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/{cell_type}_fragments.bed'
    output: 
        xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.xls",
        narrow_peak = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.narrowPeak"
    params:
        out_dir = work_dir + "/data/celltypes/{cell_type}"
    resources:
        mem_mb=200000, runtime=2880
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
        fragment_file = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
        consensus_bed = work_dir + '/data/consensus_regions.bed'
    output:
        cistopic_objects = data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
        cistopic_adata=data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample='{sample}',
        batch='{batch}',
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
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
            zip,
            sample=samples,
            batch=batches
            ),
        cistopic_adata=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad',
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
        merged_atac_anndata = work_dir + '/atlas/04_modeled_anndata_atac.h5ad',
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad'
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
        pseudo_fragment_file = work_dir + '/data/celltypes/{cell_type}/{cell_type}_fragments.bed'
    output:
        celltype_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_bigwig.bw',
        celltype_normalized_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_normalized_bigwig.bw'
    resources:
        mem_mb=1500000, runtime=960,  slurm_partition='largemem'
    singularity:
        envs['atac_fragment']
    script:
        'scripts/atac_bigwig.py'

rule celltype_bed:
    input:
        xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.xls",
        blacklist = work_dir + '/input/hg38-blacklist.bed'
    singularity:
        envs['atac_fragment']
    output:
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed'
    script:
        'scripts/MACS_to_bed.py'

rule annotate_bed:
    input:
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed'
    output:
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed'
    resources:
        runtime=30, mem_mb=50000, 
    shell:
        'module load homer;annotatePeaks.pl {input.cell_bedfile} hg38 > {output.cell_annotated_bedfile}'

rule export_atac_cell:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed',
        fragment_files=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample_key = sample_key,
        seq_batch_key = seq_batch_key,
        disease_param = disease_param,
        covariates = design_covariates,
        samples=samples,
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    threads:
        8
    resources:
        runtime=2880, mem_mb=400000, slurm_partition='largemem'
    script:
        'scripts/atac_by_celltype.py'

rule export_celltypes:
    input:
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
        separating_cluster = 'cell_type',
        design_factors = design_covariates,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/atac_DAR.py'
   
rule atac_coaccessibilty:
    input:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
    output:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac_circe.h5ad',
        circe_network = work_dir+'/data/celltypes/{cell_type}/circe_network_{cell_type}.csv'
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

rule fragments_pseudobulk_cell_disease:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        fragment_file=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        pseudo_fragment_files = expand(
            work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bed',
            cell_type=cell_types,
            disease=diseases + [control]
        )
    params:
        pseudobulk_param = 'cell_type',
        samples=samples,
        sample_param_name = sample_key,
        cell_types = cell_types,
        diseases = diseases + [control],
        disease_param = disease_param
    singularity:
        envs['atac_fragment']
    threads:
        64
    resources:
        runtime=960, mem_mb=3000000, disk_mb=500000, slurm_partition='largemem'
    script:
        'scripts/cell_disease_pseudobulk.py'

rule MACS2_peak_cell_disease:
    input:
        pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bed'
    output: 
        xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.xls",
        narrow_peak = work_dir + "/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.narrowPeak"
    params:
        out_dir = work_dir + "/data/celltypes/{cell_type}",
        cell_type = lambda wildcards: wildcards.cell_type,
        disease = lambda wildcards: wildcards.disease
    resources:
        mem_mb=200000, runtime=960
    singularity:
        envs['scenicplus']
    shell:
        "macs2 callpeak --treatment {input.pseudo_fragment_files} --name {wildcards.cell_type}_{wildcards.disease} --outdir {params.out_dir} --format BEDPE --gsize hs --qvalue 0.001 --nomodel --shift 73 --extsize 146 --keep-dup all"

rule create_bigwig_cell_disease:
    input:
        pseudo_fragment_file = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bed'
    output:
        celltype_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_bigwig.bw',
        celltype_normalized_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_normalized_bigwig.bw'
    resources:
        mem_mb=1000000, runtime=400, slurm_partition='largemem'
    singularity:
        envs['atac_fragment']
    script:
        'scripts/atac_bigwig.py'

rule celltype_bed_cell_disease:
    input:
        xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.xls",
        blacklist = work_dir + '/data/CARD_singlecell/SN_atlas/input/hg38-blacklist.bed'
    singularity:
        envs['atac_fragment']
    output:
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.bed'
    script:
        'scripts/MACS_to_bed.py'

rule annotate_bed_cell_disease:
    input:
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.bed'
    output:
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_annotated_peaks.bed'
    resources:
        runtime=30, mem_mb=50000, 
    shell:
        'module load homer;annotatePeaks.pl {input.cell_bedfile} hg38 > {output.cell_annotated_bedfile}'

rule export_atac_cell_disease:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.bed',
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_annotated_peaks.bed',
        fragment_files=expand(
            data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
            zip,
            sample=samples,
            batch=batches
            )
    output:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample_key = sample_key,
        seq_batch_key = seq_batch_key,
        disease_param = disease_param,
        covariates = design_covariates,
        samples=samples,
        cell_type = lambda wildcards: wildcards.cell_type,
        disease = lambda wildcards: wildcards.disease
    threads:
        8
    resources:
        runtime=1440, mem_mb=400000, slurm_partition='largemem'
    script:
        'scripts/atac_by_celltype.py'

rule atac_coaccessibilty_cell_disease:
    input:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_atac.h5ad'
    output:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_atac_circe.h5ad',
        circe_network = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_circe_network.csv'
    params:
        cell_type = lambda wildcards: wildcards.cell_type
    singularity:
        envs['circe']
    threads:
        8
    resources:
        runtime=600, mem_mb=400000, slurm_partition='largemem'
    script:
        'scripts/circe_by_celltype.py'

rule leiden_DGE:
    input:
        rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/significant_genes/rna/leiden/rna_leiden_{cell_type}_{disease}_DGE.csv',
        output_figure = work_dir + '/figures/leiden/rna_leiden_{cell_type}_{disease}_DGE.svg',
        celltype_pseudobulk = work_dir+'/data/celltypes/leiden/rna_leiden_{cell_type}_{disease}_pseudobulk.csv'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3],
        sample_key=sample_key,
        design_factors = design_covariates,
        separating_cluster = 'leiden_2'
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/rna_DGE.py'

rule leiden_disease_vs_disease_DGE:
    input:
        rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/significant_genes/rna/leiden/rna_{cell_type}_PD_vs_{disease}_DGE.csv',
        output_figure = work_dir + '/figures/leiden/rna_leiden_{cell_type}_PD_vs_{disease}_DGE.svg',
        celltype_pseudobulk = work_dir+'/data/celltypes/leiden/rna_leiden_{cell_type}_PD_vs_{disease}_pseudobulk.csv'
    params:
        disease_param = disease_param,
        control = 'PD',
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3],
        sample_key=sample_key,
        design_factors = design_covariates,
        separating_cluster = 'leiden_2'
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/rna_DGE.py'


rule DAR_disease_vs_disease_leiden:
    input:
        atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad'
    output:
        output_DAR_data = work_dir+'/data/significant_genes/atac/leiden/atac_{cell_type}_PD_vs_{disease}_DAR.csv',
        output_figure = work_dir+'/figures/leiden/atac_{cell_type}_PD_vs_{disease}_DAR.svg',
        cell_specific_pseudo = work_dir+'/data/celltypes/leiden/atac_leiden_{cell_type}_PD_vs_{disease}_pseudobulk.csv'
    params:
        disease_param = disease_param,
        control = 'PD',
        sample_key=sample_key,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3],
        design_factors = design_covariates,
        separating_cluster = 'leiden_2'
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/atac_DAR.py'

rule DAR_leiden:
    input:
        atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad'
    output:
        output_DAR_data = work_dir+'/data/significant_genes/atac/leiden/atac_{cell_type}_{disease}_DAR.csv',
        output_figure = work_dir+'/figures/leiden/atac_{cell_type}_{disease}_DAR.svg',
        cell_specific_pseudo = work_dir+'/data/celltypes/leiden/atac_leiden_{cell_type}_{disease}_pseudobulk.csv'
    params:
        disease_param = disease_param,
        control = control,
        sample_key=sample_key,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3],
        design_factors = [],
        separating_cluster = 'leiden_2'
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/atac_DAR.py'

rule motif_enrichment:
    input:
        atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad',
        ref_genome = reference_genome,
        TF_motifs = work_dir + '/input/jaspar_2024_hsapiens.meme'
    output:
        motif_enrichment = work_dir+'/data/motif_enrichment.csv'
    params:
        control = control,
        cell_type = 'cell_type',
        disease_param = disease_param
    singularity:
        envs['snapatac2']
    resources:
        runtime=240, disk_mb=300000, mem_mb=200000
    script:
        'scripts/atac_motif_enrichment.py'

rule differential_motif_enrichment:
    input:
        output_DAR_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
        cell_type_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad',
        TF_motifs = work_dir + '/input/jaspar_2024_hsapiens.meme',
        ref_genome = reference_genome
    output:
        differential_motif_dataframe = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_differential_motif.csv'
    singularity:
        envs['snapatac2']
    resources:
        runtime=240, disk_mb=300000, mem_mb=200000
    script:
        'scripts/differential_motif_enrichment.py'

rule DAR_CCAN_modules:
    input:
        celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac_circe.h5ad',
        output_DAR_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv'
    output:
        output_DAR_CCAN_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_CCAN_DAR.csv'
    singularity:
        envs['circe']
    resources:
        runtime=240, disk_mb=300000, mem_mb=200000
    script:
        'scripts/atac_DAR_CCANs.py'

rule disease_gsea:
    input:
        adata_path =  work_dir+'/data/celltypes/{cell_type}/rna.h5ad',
        ontologies = work_dir+'/input/ontologies.csv'
    output:
        cell_disease_GSEA =  work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_GSEA_genes.csv'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['great_gsea']
    threads:
        64
    resources:
        runtime=960, mem_mb=1000000, slurm_partition='largemem' 
    script:
        'scripts/rna_GSEA.py'

rule disease_great:
    input:
        DAR_path =  work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
        tss_file =  work_dir+'/input/tss_from_great.bed',
        chr_sizes_file =  work_dir+'/input/chr_size.bed',
        annotation_file =  work_dir+'/input/ontologies.csv',
    output:
        cell_disease_peaks = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_DAR_peaks.bed',
        cell_disease_GREAT = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_GREAT_peaks.csv'
    singularity:
        envs['great_gsea']
    resources:
        runtime=1440
    script:
        'scripts/atac_GREAT.py'

rule celltype_disease_bed2bam:
    input:
        bed = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bed',
        ref_genome_length = genome_length
    output:
        bam = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bam'
    singularity:
        envs['atac_fragment']
    resources:
        runtime=2880, mem_mb=300000
    shell:
        'bedToBam -i {input.bed} -g {input.ref_genome_length} > {output.bam}'

rule celltype_disease_ATACorrect:
    input:
        bam = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bam',
        blacklist = work_dir + '/input/hg38-blacklist.bed',
        cell_type_peaks = work_dir+'/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        ref_genome = reference_genome
    output:
        corrected_bigwig = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_corrected.bw'
    params:
        ATACorrect_outdir = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect',
        prefix = '{cell_type}_{disease}'
    singularity:
        envs['atac_fragment']
    threads:
        64
    resources:
        runtime=960, mem_mb=300000
    shell:
        'TOBIAS ATACorrect --bam {input.bam} --genome {input.ref_genome} --blacklist {input.blacklist} --peaks {input.cell_type_peaks} --outdir {params.ATACorrect_outdir} --prefix {params.prefix} --cores {threads}'

rule celltype_disease_score_bigwig:
    input:
        corrected_bigwig = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_corrected.bw',
        regions = work_dir+'/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
    output:
        footprinted_bigwig = work_dir+'/data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_footprints.bw'
    singularity:
        envs['atac_fragment']
    threads:
        64
    resources:
        runtime=960, mem_mb=300000
    shell:
        'TOBIAS FootprintScores --signal {input.corrected_bigwig} --regions {input.regions} --output {output.footprinted_bigwig} --cores {threads}'

rule control_comparison_score_bigwig:
    input:
        corrected_bigwig = work_dir+'/data/celltypes/{cell_type}/{cell_type}_control_ATACorrect/{cell_type}_control_corrected.bw',
        regions = work_dir+'/data/consensus_regions.bed'
    output:
        control_footprint_bigwig = work_dir+'/data/celltypes/{cell_type}/{cell_type}_control_ATACorrect/{cell_type}_control_comparison_footprints.bw'
    singularity:
        envs['atac_fragment']
    threads:
        64
    resources:
        runtime=960, mem_mb=300000
    shell:
        'TOBIAS FootprintScores --signal {input.corrected_bigwig} --regions {input.regions} --output {output.control_footprint_bigwig} --cores {threads}'