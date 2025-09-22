#IMPORTS
import os
import pandas as pd

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

"""File locations"""
data_dir            = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/'      # Define the data directory, explicitly
# work_dir            = os.getcwd()                                           # Define the working directory, explictly as the directory of this pipeline
work_dir            = '/vf/users/CARD_singlecell/SN_control_atlas/scMAVERICS'
metadata_table      = work_dir + '/input/SN_control_samples-EF_fixed.csv'   # Define where the metadata data exists for each sample to be processed
gene_markers_file   = work_dir + '/input/celltype_markers_dict_reduced.csv' # Define where celltypes/cell marker gene 
# # TEMP
# my_batch_list = ['2', '3', '2', '2', '3', '2', '3', '2', '0', '2', '6']
# my_sample_list = ['2539', '5931', '1424', '1297', '2781', '1745', '4276', '2420', '1560', '1615', 'PT13122']
# # </TEMP>

"""Metadata parameters"""
seq_batch_key = 'Use_batch'                                     # Key for sequencing batch, used for directory search
sample_key = 'Sample'                                           # Key for samples, required in aggregating while preserving sample info
batches = pd.read_csv(metadata_table)[seq_batch_key].tolist()   # Read in the list of batches and samples
samples = pd.read_csv(metadata_table)[sample_key].tolist()
disease_param = 'age_bin'                                       # Name of the disease parameter
control = 'young'                                               # Define disease states
diseases = ['old']                                              # Disease states to compare, keep as list of strings, unnecessary 
cell_types = pd.read_csv(gene_markers_file)['cell type']        # Define the cell types to look for, from gene marker file
# design_covariates = [seq_batch_key, 'Age', 'Sex']               # Design factors/covariates for DGEs and DARs
design_covariates = ['Sex']                                     # Design factors/covariates for DGEs and DARs

"""Quality control thresholds"""
mito_percent_thresh = 15    # Maximum percent of genes in a cell that can be mitochondrial
ribo_percent_thresh = 10    # Maximum percent of genes in a cell that can be ribosomal
doublet_thresh      = 0.15  # Maximum doublet score for a cell, computed by scrublet
min_genes_per_cell  = 250   # Minimum number of unique genes in a cell
min_peak_counts     = 500   # Minimum number of fragments per cell


"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'snapatac2'     : 'envs/snapatac2.sif',
    'singlecell'    : 'envs/single_cell_gpu.sif',
    'scenicplus'    : 'envs/scenicplus.sif',
    'decoupler'     : 'envs/decoupler.sif',
    'circe'         : 'envs/circe.sif',
    'atac_fragment' : 'envs/atac_fragment.sif'
    }

rule all:
    input:
        # # EF TEST
        # atac_anndata = expand(
        #     # data_dir + '{sample}/01_{sample}_anndata_object_atac.h5ad',
        #     data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad',
        #     zip,
        #     sample = samples,
        #     sample_key = sample_key,
        #     batch = batches,
        #     ),
        # atac_anndata = expand(
        #     # data_dir + '{sample}/03_{sample}_anndata_object_atac.h5ad',
        #     data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad',
        #     zip,
        #     sample = samples,
        #     sample_key = sample_key,
        #     batch = batches
        #     ),
        # merged_rna_anndata = work_dir + '/atlas/03_filtered_anndata_rna.h5ad',
        # hvg_rna_anndata = work_dir + '/atlas/03_hvg_anndata_rna.h5ad',
        # model_history = work_dir + '/data/model_elbo/rna_model_history.csv',
        # cell_bedfile = expand(
        #     work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        #     zip,
        #     cell_type = cell_types
        #     ),
        # cell_annotated_bedfile = expand(
        #     work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed',
        #     zip,
        #     cell_type = cell_types
        #     )
        # celltype_atac = expand(
        #     work_dir + '/data/celltypes/{cell_type}/atac.h5ad',
        #     zip,
        #     cell_type = cell_types
        # ),
        # # EF - TEST - fix
        # cistopic_object = expand(
        #     data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
        #     zip,
        #     batch = my_batch_list,
        #     sample = my_sample_list
        # ),
        # cistopic_adata = expand(
        #     data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad',
        #     zip,
        #     batch = my_batch_list,
        #     sample = my_sample_list
        # ),
        # cistopic_object = expand(
        #     data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
        #     zip,
        #     batch = batches,
        #     sample = samples
        # ),
        # cistopic_adata = expand(
        #     data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad',
        #     zip,
        #     batch = batches,
        #     sample = samples
        # ),
        # # TEST - rerun just Sample 2587 - worked!
        # cistopic_object = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/batch1/Multiome/2587-ARC/outs/04_2587_cistopic_obj.pkl',
        # cistopic_adata  = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/batch1/Multiome/2587-ARC/outs/04_2587_anndata_peaks_atac.h5ad',
        # # TODO
        merged_cistopic_object = work_dir + '/data/merged_cistopic_object.pkl',
        merged_atac_anndata = work_dir + '/atlas/03_merged_cistopic_atac.h5ad'

# merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad'
# # This is the last step of the pipeline, run all the way through with this input or swap out for an intermediary file below for checkpoints

# # Uncomment to view QC data
# """genes_by_counts = work_dir+'figures/QC_genes_by_counts.png'"""
# # Uncomment when you have verified QC metrics
# rna_anndata=expand(
#             data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad', 
#             zip,
#             batch = batches,
#             sample = samples
#             ),

# # Uncomment when you want to model rna data
# merged_rna_anndata = work_dir + '/atlas/05_annotated_anndata_rna.h5ad',

# # Uncomment when you want to run DGE/DAR analysis
# output_DGE_data = expand(
#     work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
#     cell_type = cell_types,
#     disease = diseases
#     ),

# # Uncomment when you want to model ATAC-peak data
# merged_cistopic_adata = work_dir + '/atlas/05_annotated_anndata_atac.h5ad',
# merged_multiome = work_dir + '/atlas/multiome_atlas.h5mu',
# output_DAR_data = expand(
#     work_dir + '/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
#     cell_type = cell_types,
#     disease = diseases
#     ) # This is the last step of the pipeline, run all the way through with this input or swap out for an intermediary file below for checkpoints

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # #! NOTE - not sure where the cwd folder should be
# # #! WARNING - run once, then comment out!
# rule cellbender:
#     input:
#         # rna_anndata = data_dir+'{sample}/raw_feature_bc_matrix.h5',
#         # cwd = data_dir + '{sample}/'
#         rna_anndata = data_dir + '{batch}/Multiome/{sample}-ARC/outs/raw_feature_bc_matrix.h5',
#         cwd = data_dir + '{batch}/Multiome/{sample}-ARC/'
#     output:
#         # rna_anndata = data_dir + '{sample}/cellbender_gex_counts_filtered.h5'
#         rna_anndata = data_dir + '{batch}/Multiome/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
#     params:
#         sample = '{sample}',
#         batch = '{batch}'
#     resources:
#         runtime=2880, mem_mb=300000, gpu=1, gpu_model='v100x'
#     shell:
#         work_dir + '/scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rule rna_preprocess:
#     input:
#         metadata_table = metadata_table,
#         # rna_anndata = data_dir + '{sample}/cellbender_gex_counts_filtered.h5'
#         rna_anndata = data_dir + '{batch}/Multiome/{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
#     output:
#         rna_anndata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad'
#         # rna_anndata = data_dir + '{sample}/01_{sample}_anndata_object_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     params:
#         batch = '{batch}',
#         sample = '{sample}',
#         sample_key = sample_key,
#     resources:
#         runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
#     script:
#         work_dir + '/scripts/rna_preprocess.py'

# rule merge_unfiltered:
#     input:
#         rna_anndata = expand(
#             # data_dir + '{sample}/01_{sample}_anndata_object_rna.h5ad', 
#             # data_dir + '{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad', 
#             data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_rna.h5ad', 
#             zip,
#             batch = batches,
#             sample = samples
#         )
#     output:
#         merged_rna_anndata = work_dir + '/atlas/01_merged_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     params:
#         samples = samples
#     resources:
#         runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
#     script:
#         work_dir + '/scripts/merge_anndata.py'

# rule plot_qc_rna:
#     input:
#         merged_rna_anndata  = work_dir + '/atlas/01_merged_anndata_rna.h5ad'
#     output:
#         mito_figure         = work_dir + '/figures/QC_mito_pct.png',
#         ribo_figure         = work_dir + '/figures/QC_ribo_pct.png',
#         gene_counts_figure  = work_dir + '/figures/QC_gene_counts.png',
#         doublet_figure      = work_dir + '/figures/QC_doublet.png',
#         genes_by_counts     = work_dir + '/figures/QC_genes_by_counts.png'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
#     params:
#         mito_percent_thresh = mito_percent_thresh,
#         doublet_thresh      = doublet_thresh,
#         min_genes_per_cell  = min_genes_per_cell,
#         ribo_percent_thresh = ribo_percent_thresh,
#         sample_key          = sample_key,
        
#     script:
#         work_dir + '/scripts/plot_qc_metrics.py'

# rule filter_rna:
#     input:        
#         rna_anndata = data_dir + '{sample}/01_{sample}_anndata_object_rna.h5ad'
#     output:
#         rna_anndata = data_dir + '{sample}/02_{sample}_anndata_filtered_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     params:
#         mito_percent_thresh     = mito_percent_thresh,
#         doublet_thresh          = doublet_thresh,
#         min_genes_per_cell      = min_genes_per_cell,
#         ribo_percent_thresh     = ribo_percent_thresh
#     resources:
#         runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
#     script: 
#         work_dir + '/scripts/rna_filter.py'

# rule merge_filtered_rna:
#     input:
#         rna_anndata = expand(
#             data_dir + '{sample}/02_{sample}_anndata_filtered_rna.h5ad', 
#             zip,
#             batch = batches,
#             sample = samples
#             )
#     output:
#         merged_rna_anndata = work_dir + '/atlas/02_filtered_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     params:
#         samples = samples
#     resources:
#         runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
#     script:
#         work_dir + '/scripts/merge_anndata.py'

# rule atac_preprocess:
#     input:
#         # fragment_file = data_dir + '{sample}/atac_fragments.tsv.gz'
#         # fragment_file = data_dir + 'batch{batch}/Multiome/{sample}-ARC/atac_fragments.tsv.gz'
#         fragment_file = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz'
#     output:
#         # atac_anndata = data_dir + '{sample}/01_{sample}_anndata_object_atac.h5ad'
#         atac_anndata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad'
#     singularity:
#         envs['snapatac2']
#     resources:
#         runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
#     script:
#         work_dir+'/scripts/atac_preprocess.py'

# #! NOTE - don't need this rule
# rule merge_unfiltered_atac:
#     input:
#         metadata_table = metadata_table,
#         atac_anndata = expand(
#             # data_dir+'{sample}/01_{sample}_anndata_object_atac.h5ad',
#             data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad',
#             zip,
#             batch=batches,
#             sample=samples
#             ),
#     output:
#         merged_atac_anndata = work_dir+'/atlas/01_merged_anndata_atac.h5ad'
#     singularity:
#         envs['snapatac2']
#     params:
#         samples=samples,
#         sample_key=sample_key
#     resources:
#         runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
#     script:
#         work_dir+'/scripts/merge_atac.py'

# #* NOTE - don't need to run, should just plot them manually
# rule plot_qc_atac:
#     input:
#         atac_anndata=data_dir+'{sample}/01_{sample}_anndata_object_rna.h5ad' #? Q) shouldn't this be _atac.h5ad?
#     singularity:
#         envs['snapatac2']
#     resources:
#         runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/atac_plot_qc.py'

# rule filter_atac:
#     input:
#         # rna_anndata = data_dir+'{sample}/02_{sample}_anndata_filtered_rna.h5ad',
#         # atac_anndata = data_dir + '{sample}/01_{sample}_anndata_object_atac.h5ad',
#         rna_anndata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/02_{sample}_anndata_filtered_rna.h5ad',
#         atac_anndata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad'
#     output:
#         # atac_anndata = data_dir + '{sample}/03_{sample}_anndata_object_atac.h5ad',
#         # rna_anndata = data_dir + '{sample}/03_{sample}_anndata_filtered_rna.h5ad'
#         atac_anndata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_object_atac.h5ad',
#         rna_anndata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad'
#     singularity:
#         envs['snapatac2']
#     resources:
#         runtime=30, mem_mb=50000, slurm_partition='quick'
#     script:
#         work_dir+'/scripts/atac_filter.py'

# rule merge_multiome_rna:
#     input:
#         rna_anndata=expand(
#             data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/03_{sample}_anndata_filtered_rna.h5ad', 
#             zip,
#             batch=batches,
#             sample=samples
#             )
#     output:
#         merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     params:
#         samples=samples
#     resources:
#         runtime=120, mem_mb=300000, disk_mb=10000, slurm_partition='norm'
#         #slurm_partition='largemem' 
#     script:
#         work_dir+'/scripts/merge_anndata.py'

# rule feature_selection:
#     input:
#         merged_rna_anndata = work_dir+'/atlas/03_filtered_anndata_rna.h5ad'
#     output:
#         hvg_rna_anndata = work_dir+'/atlas/03_hvg_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=360, mem_mb=1500000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/feature_selection.py'

# rule rna_model:
#     input:
#         hvg_rna_anndata = work_dir + '/atlas/03_hvg_anndata_rna.h5ad'
#     output:
#         hvg_rna_anndata = work_dir + '/atlas/04_modeled_hvg_anndata_rna.h5ad',
#         model_history = work_dir + '/data/model_elbo/rna_model_history.csv'
#     params:
#         model = work_dir + '/data/models/rna/',
#         sample_key = sample_key
#     threads:
#         64
#     resources:
#         runtime=2880, mem_mb=300000, gpu=4, gpu_model='v100x'
#     shell:
#         'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model}'

# rule UMAP:
#     input:
#         merged_rna_anndata = work_dir + '/atlas/03_filtered_anndata_rna.h5ad',
#         hvg_rna_anndata = work_dir + '/atlas/04_modeled_hvg_anndata_rna.h5ad'
#     output:
#         merged_rna_anndata = work_dir + '/atlas/04_modeled_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=1440, mem_mb=1000000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/scVI_to_UMAP.py'

# rule first_pass_annotate:
#     input:
#         merged_rna_anndata = work_dir+'/atlas/04_modeled_anndata_rna.h5ad',
#         # gene_markers = work_dir+'/input/first_pass_genes.csv'
#         gene_markers = work_dir+'/input/celltype_markers_dict_reduced.csv'
#     output:
#         merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad',
#         cell_annotate = work_dir+'/data/first_pass_genes.csv'
#     params:
#         seq_batch_key = seq_batch_key
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=480, mem_mb=1500000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/annotate.py'

# rule cluster_based_QC:
#     input:
#         merged_rna_anndata = work_dir + '/atlas/05_annotated_anndata_rna.h5ad'
#         # merged_rna_anndata = work_dir + '/atlas/05_annotated_anndata_rna-EF_manual.h5ad'
#     output:
#         merged_rna_anndata = work_dir + '/atlas/05_QC_filtered_anndata_rna.h5ad',
#         course_celltype = work_dir + '/figures/first_pass_RNA_UMAP_celltype.svg',
#         course_counts = work_dir + '/figures/first_pass_RNA_num_genes_celltype.svg'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=240, mem_mb=1500000, slurm_partition='largemem'
#     script:
#         work_dir + '/scripts/cluster_based_QC.py'

# rule filtered_feature_selection:
#     input:
#         merged_rna_anndata = work_dir+'/atlas/05_QC_filtered_anndata_rna.h5ad'
#     output:
#         hvg_rna_anndata = work_dir+'/atlas/05_hvg_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=360, mem_mb=1500000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/feature_selection.py'

# rule rna_polish_model:
#     input:
#         hvg_rna_anndata = work_dir+'/atlas/05_hvg_anndata_rna.h5ad'
#     output:
#         hvg_rna_anndata = work_dir+'/atlas/05_modeled_hvg_anndata_rna.h5ad',
#         model_history = work_dir+'/data/model_elbo/rna_model_v2_history.csv'
#     params:
#         model = work_dir+'/data/models/rna_polish/',
#         sample_key = sample_key
#     threads:
#         64
#     resources:
#         runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
#     shell:
#         'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model}'

# rule filtered_UMAP:
#     input:
#         merged_rna_anndata = work_dir + '/atlas/05_QC_filtered_anndata_rna.h5ad',
#         hvg_rna_anndata = work_dir + '/atlas/05_modeled_hvg_anndata_rna.h5ad'
#     output:
#         merged_rna_anndata = work_dir + '/atlas/06_polished_anndata_rna.h5ad'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=1440, mem_mb=1000000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/scVI_to_UMAP.py'

# rule second_pass_annotate:
#     input:
#         merged_rna_anndata = work_dir+'/atlas/06_polished_anndata_rna.h5ad',
#         gene_markers = gene_markers_file
#     output:
#         merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
#         cell_annotate = work_dir+'/data/rna_cell_annot.csv'
#     params:
#         seq_batch_key = seq_batch_key
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=240, mem_mb=1500000, slurm_partition='largemem'
#     script:
#         work_dir+'/scripts/annotate.py'

# rule DGE:
#     input:
#         # rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad'
#         rna_anndata = work_dir + '/atlas/07_polished_anndata_rna_recalc_missing_rows.h5ad'
#     output:
#         output_DGE_data = work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
#         output_figure = work_dir + '/figures/{cell_type}/rna_{cell_type}_{disease}_DGE.svg',
#         celltype_pseudobulk = work_dir+'/data/celltypes/{cell_type}/rna_{cell_type}_{disease}_pseudobulk.csv'
#     params:
#         disease_param = disease_param,
#         control = control,
#         disease = lambda wildcards, output: output[0].split("_")[-2],
#         cell_type = lambda wildcards, output: output[0].split("_")[-3],
#         sample_key=sample_key,
#         design_factors = design_covariates
#     singularity:
#         envs['decoupler']
#     threads:
#         64
#     resources:
#         runtime=1440, disk_mb=200000, mem_mb=200000
#     script:
#         'scripts/rna_DGE.py'

#! WARNING - from here on, replace the 07_polished_anndata_rna.h5ad with 07_polished_anndata_rna_recalc_missing_rows.h5ad
# rule cistopic_pseudobulk:
#     input:
#         # merged_rna_anndata = work_dir + '/atlas/07_annotated_anndata_rna.h5ad',
#         # merged_rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad',
#         merged_rna_anndata = work_dir + '/atlas/07_polished_anndata_rna_recalc_missing_rows.h5ad',
#         fragment_file = expand(
#             data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
#             zip,
#             sample = samples,
#             batch = batches
#             )
#     output:
#         pseudo_fragment_files = expand(
#             work_dir + '/data/celltypes/{cell_types}/{cell_types}_fragments.bed',
#             cell_types=cell_types)
#     params:
#         pseudobulk_param = 'cell_type',
#         samples = samples,
#         sample_param_name = sample_key,
#         cell_types = cell_types
#     singularity:
#         envs['atac_fragment']
#     threads:
#         64
#     resources:
#         runtime=960, mem_mb=3000000, disk_mb=500000, slurm_partition='largemem'
#     script:
#         'scripts/cistopic_pseudobulk.py'

# rule MACS2_peak_call:
#     input:
#         pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/{cell_type}_fragments.bed'
#     output: 
#         xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.xls",
#         narrow_peak = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.narrowPeak"
#     params:
#         out_dir = work_dir + "/data/celltypes/{cell_type}"
#     resources:
#         mem_mb=200000, runtime=960
#     singularity:
#         envs['scenicplus']
#     shell:
#         "macs2 callpeak --treatment {input.pseudo_fragment_files} --name {wildcards.cell_type} --outdir {params.out_dir} --format BEDPE --gsize hs --qvalue 0.001 --nomodel --shift 73 --extsize 146 --keep-dup all"

# rule consensus_peaks:
#     input:
#         narrow_peaks = expand(
#             work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.narrowPeak",
#             cell_type = cell_types
#             )
#     output:
#         consensus_bed = work_dir + '/data/consensus_regions.bed'
#     singularity:
#         envs['scenicplus']
#     resources:
#         runtime=960, mem_mb=100000
#     script:
#         'scripts/MACS_consensus.py'

# rule cistopic_create_objects:
#     input:
#         # merged_rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad',
#         merged_rna_anndata = work_dir + '/atlas/07_polished_anndata_rna_recalc_missing_rows.h5ad',
#         # fragment_file = data_dir + '{sample}/atac_fragments.tsv.gz',
#         fragment_file = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
#         consensus_bed = work_dir + '/data/consensus_regions.bed'
#     output:
#         # cistopic_object = data_dir + '{sample}/04_{sample}_cistopic_obj.pkl',
#         # cistopic_adata = data_dir + '{sample}/04_{sample}_anndata_peaks_atac.h5ad'
#         cistopic_object = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
#         cistopic_adata = data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad'
#     singularity:
#         envs['scenicplus']
#     params:
#         batch = '{batch}',
#         sample = '{sample}',
#         sample_key = sample_key,
#         seq_batch_key = seq_batch_key,
#         disease_param = disease_param
#     resources:
#         # runtime=120, mem_mb=250000, slurm_partition='quick'
#         # runtime=120, mem_mb=250000, slurm_partition='norm'
#         runtime=120, mem_mb=240000, slurm_partition='quick'
#     threads:
#         16
#     script:
#         'scripts/cistopic_create_object.py'

rule cistopic_merge_objects:
    input:
        # merged_rna_anndata = work_dir+'/atlas/05_polished_anndata_rna.h5ad',
        merged_rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad',
        cistopic_objects = expand(
            # data_dir+'{sample}/04_{sample}_cistopic_obj.pkl',
            data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_cistopic_obj.pkl',
            zip,
            sample = samples,
            batch = batches
            ),
        rna_anndata = expand(
            # data_dir+'{sample}/04_{sample}_anndata_peaks_atac.h5ad', 
            data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/04_{sample}_anndata_peaks_atac.h5ad',
            zip,
            sample = samples,
            batch = batches
            )
    output:
        merged_cistopic_object = work_dir + '/data/merged_cistopic_object.pkl',
        merged_atac_anndata = work_dir + '/atlas/03_merged_cistopic_atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample_key = sample_key,
        disease_param = disease_param
    resources:
        runtime=1440, mem_mb=2000000, slurm_partition='largemem'
    script:
        'scripts/merge_cistopic_and_adata.py'

# rule atac_peaks_model:
#     input:
#         merged_atac_anndata = work_dir+'/atlas/03_merged_cistopic_atac.h5ad'
#     output:
#         merged_atac_anndata = work_dir+'/atlas/04_modeled_anndata_atac.h5ad',
#         atac_model_history = work_dir+'/data/model_elbo/atac_model_history.csv'
#     params:
#         atac_model = work_dir+'/data/models/atac/',
#         sample_key = sample_key
#     threads:
#         64
#     resources:
#         runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
#     shell:
#         'scripts/atac_model.sh {input.merged_atac_anndata} {params.sample_key} {output.atac_model_history} {output.merged_atac_anndata} {params.atac_model}'

# rule multiome_output:
#     input:
#         merged_atac_anndata = work_dir + '/atlas/04_modeled_anndata_atac.h5ad',
#         merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad'
#     output:
#         merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
#     singularity:
#         envs['singlecell']
#     resources:
#         runtime=120, mem_mb=300000, slurm_partition='quick' 
#     script:
#         'scripts/merge_muon.py'

# rule create_bigwig:
#     input:
#         pseudo_fragment_file = work_dir + '/data/celltypes/{cell_type}/{cell_types}_fragments.bed'
#     output:
#         celltype_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_bigwig.bw',
#         celltype_normalized_bigwig = work_dir + '/data/celltypes/{cell_type}/{cell_type}_normalized_bigwig.bw'
#     resources:
#         mem_mb=300000, runtime=180, slurm_partition='quick'
#     singularity:
#         envs['atac_fragment']
#     script:
#         'scripts/atac_bigwig.py'

# rule celltype_bed:
#     input:
#         xls = work_dir + "/data/celltypes/{cell_type}/{cell_type}_peaks.xls",
#     singularity:
#         envs['atac_fragment']
#     output:
#         cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed'
#     script:
#         'scripts/MACS_to_bed.py'

# rule annotate_bed:
#     input:
#         cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed'
#     output:
#         cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed'
#     resources:
#         runtime=30, mem_mb=50000, 
#     shell:
#         'module load homer;annotatePeaks.pl {input.cell_bedfile} hg38 > {output.cell_annotated_bedfile}'

rule export_atac_cell:
    input:
        merged_rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad',
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed',
        fragment_files = expand(
            data_dir + 'batch{batch}/Multiome/{sample}-ARC/outs/atac_fragments.tsv.gz',
            zip,
            sample = samples,
            batch = batches
            )
    output:
        celltype_atac = work_dir + '/data/celltypes/{cell_type}/atac.h5ad'
    singularity:
        envs['scenicplus']
    params:
        sample_key = sample_key,
        seq_batch_key = seq_batch_key,
        disease_param = disease_param,
        samples = samples,
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    threads:
        8
    resources:
        runtime=2880, mem_mb=400000, slurm_partition='largemem'
    script:
        'scripts/atac_by_celltype.py'

# """rule export_celltypes:
#     input:
#         merged_multiome = work_dir+'/atlas/multiome_atlas.h5mu'
#     output:
#         celltype_atac = work_dir+'/data/celltypes/{cell_type}/consensus_peaks_atac.h5ad',
#         celltype_rna = work_dir+'/data/celltypes/{cell_type}/rna.h5ad'
#     params:
#         cell_type = lambda wildcards, output: output[0].split('/')[-2]
#     singularity:
#         envs['singlecell']
#     threads:
#         8
#     resources:
#         runtime=120, mem_mb=300000
#     script:
#         'scripts/export_celltype.py'
# """

# rule DAR:
#     input:
#         atac_anndata = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
#     output:
#         output_DAR_data = work_dir+'/data/significant_genes/atac/atac_{cell_type}_{disease}_DAR.csv',
#         output_figure = work_dir+'/figures/{cell_type}/atac_{cell_type}_{disease}_DAR.svg',
#         cell_specific_pseudo = work_dir+'/data/celltypes/{cell_type}/atac_{disease}_pseudobulk.csv'
#     params:
#         disease_param = disease_param,
#         control = control,
#         disease = lambda wildcards, output: output[0].split("_")[-2],
#         cell_type = lambda wildcards, output: output[0].split("_")[-3],
#         design_factors = design_covariates
#     singularity:
#         envs['decoupler']
#     threads:
#         64
#     resources:
#         runtime=1440, disk_mb=200000, mem_mb=200000
#     script:
#         'scripts/atac_DAR.py'
   
# rule atac_coaccessibilty:
#     input:
#         celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac.h5ad'
#     output:
#         celltype_atac = work_dir+'/data/celltypes/{cell_type}/atac_circe.h5ad',
#         circe_network = work_dir+'/data/celltypes/{cell_type}/circe_network_{cell_type}.csv'
#     params:
#         cell_type = lambda wildcards, output: output[0].split('/')[-2]
#     singularity:
#         envs['circe']
#     threads:
#         8
#     resources:
#         runtime=2880, mem_mb=1500000, slurm_partition='largemem'
#     script:
#         'scripts/circe_by_celltype.py'
