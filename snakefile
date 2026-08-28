import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

configfile: "config.yaml"
"""File locations"""
data_dir = config['data_dir'] # Define the data directory, explicitly
work_dir = config['work_dir'] # Define the working directory, explictly as the directory of this pipeline
metadata_table = work_dir+config['metadata'] # Define where the metadata data exists for each sample to be processed
gene_markers_file = work_dir+config['gene_list'] # Define where celltypes/cell marker gene 
cell_cycle_gene_file = work_dir+config['cell_cycle_genes']
gene_info = work_dir+config['gene_info']
gene_tss = work_dir+config['gene_tss']
motifs = work_dir + config['motif_file']
cell_cycle_gene_file = work_dir + config['lab_cell_cycle_genes']

"""Metadata parameters"""
seq_batch_key = config['seq_batch_key'] # Key for sequencing batch, used for directory search
sample_key = config['sample_key'] # Key for samples, required in aggregating while preserving sample info

batches = pd.read_csv(metadata_table)[seq_batch_key].tolist() # Read in the list of batches and samples
samples = pd.read_csv(metadata_table)[sample_key].tolist()

disease_param = config['disease_param'] # Name of the disease parameter
control = config['control_key'] # Define disease states
diagnoses = config['diagnoses'] # Disease states to compare, keep as list of strings, unnecessary 
#disease_comparisons = ['control vs. PD', 'control vs. DLB', 'PD vs. DLB']

cell_types = pd.read_csv(gene_markers_file)['cell type'] # Define the cell types to look for, from gene marker file
design_covariates = config['covariates'] # Design factors/covariates for DGEs and DARs
reference_genome = config['reference_genome']
genome_length = config['genome_length']

"""Quality control thresholds"""
mito_percent_thresh = config['mito_thresh']# Maximum percent of genes in a cell that can be mitochondrial
ribo_percent_thresh = config['ribo_thresh'] # Maximum percent of genes in a cell that can be ribosomal
doublet_thresh = config['doublet_thresh'] # Maximum doublet score for a cell, computed by scrublet
min_genes_per_cell = config['min_genes'] # Minimum number of unique genes in a cell
min_peak_counts = config['min_peaks'] # Minimum number of fragments per cell
min_tsse = config['min_tsse'] # Minimum transcription start site enrichment

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'singlecell': 'envs/single_cell_gpu.sif',
    'tobias': 'envs/tobias.sif',
    'dreampy': 'envs/dreampy.sif',
    'multiome': 'envs/multiome.sif',
    'scenic': 'envs/scenicplus.sif'
    }

rule all:
    input:
        merged_rna_anndata = work_dir+'atlas/05_QC_filtered_anndata_rna.h5ad',

"""========================================================================="""
"""                                RNA portion                              """
"""========================================================================="""

rule cellbender:
    input:
        rna_anndata =data_dir+'{sample}/outs/raw_feature_bc_matrix.h5',
        cwd = data_dir+'{sample}-ARC/outs/'
    output:
        rna_anndata = data_dir+'{sample}-ARC/outs/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=1440, mem_mb=200000, gpu=1, gpu_model='v100x'
    shell:
        work_dir+'scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'

rule rna_preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'{sample}/outs/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'{sample}/outs/01_{sample}_anndata_object_rna.h5ad'
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
        rna_anndata = data_dir+'{sample}/outs/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'{sample}/outs/02_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['multiome']
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
    script: 
        work_dir+'scripts/rna_filter.py'

rule rna_merge:
    input:
        rna_anndata=expand(
            data_dir+'{sample}/outs/02_{sample}_anndata_filtered_rna.h5ad', 
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
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/rna_model.sh {input.hvg_rna_anndata} {params.sample_key} {output.model_history} {output.hvg_rna_anndata} {params.model} {params.random_number_seed} {params.num_layers} {params.num_latent} {params.max_epoch} {params.machine_type}'


rule rna_latent_transfer:
    input:
        merged_rna_anndata = work_dir + '/atlas/02_filtered_anndata_rna.h5ad',
        hvg_rna_anndata = work_dir + '/atlas/04_modeled_hvg_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir + '/atlas/04_modeled_anndata_rna.h5ad'
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
        course_celltype = work_dir + '/figures/first_pass_RNA_UMAP_celltype.svg',
        course_counts = work_dir + '/figures/first_pass_RNA_num_genes_celltype.svg'
    singularity:
        envs['multiome']
    resources:
        runtime=240, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir + '/scripts/rna_cluster_based_QC.py'


"""========================================================================="""
"""                               ATAC portion                              """
"""========================================================================="""

rule atac_preprocess:
    input:
        fragment_file=data_dir+'{sample}/outs/atac_fragments.tsv.gz'
    output:
        atac_anndata=data_dir+'{sample}/outs/01_{sample}_anndata_object_atac.h5ad'
    singularity:
        envs['multiome']
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'scripts/atac_preprocess.py'

rule atac_filter:
    input:
        atac_anndata = data_dir+'{sample}/outs/01_{sample}_anndata_object_atac.h5ad'
    output:
        atac_anndata = data_dir+'{sample}/outs/02_{sample}_anndata_filtered_atac.h5ad'
    singularity:
        envs['multiome']
    params:
        min_peak_counts = min_peak_counts,
        min_tsse = min_tsse
    script:
        work_dir+'scripts/atac_filter.py'

rule atac_merge:
    input:
        fragments=expand(
            data_dir+'{sample}/outs/atac_fragments.tsv.gz', 
            zip,
            batch=working_batches,
            sample=working_samples
            )
    output:
        atac_anndata = expand(
            data_dir+'{sample}/outs/02_{sample}_anndata_filtered_atac.h5ad',
            zip,
            batch=working_batches,
            sample=working_samples
            ),
        temp_file = work_dir+'/atlas/temp_filtered_anndata_atac.h5ad',
        merged_atac_anndata = work_dir+'/atlas/02_concat_atac.h5ad'
    singularity:
        envs['multiome']
    params:
        samples=working_samples
    threads:
        64
    resources:
        runtime=720, mem_mb=3000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/atac_merge.py'

rule atac_fragment_pseudobulk:
    input:
        merged_rna_anndata = work_dir+'/atlas/05_QC_filtered_anndata_rna.h5ad',
        merged_atac_anndata = work_dir + '/atlas/02_filtered_anndata_atac_backup.h5ad',
        fragment_file=expand(
            data_dir+'{batch}/Multiome/{sample}/outs/atac_fragments.tsv.gz',
            zip,
            batch=batches,
            sample=samples
            )
    output:
        pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/{cell_type}_fragments.bed'
    params:
        pseudobulk_param = 'celltype',
        samples=samples,
        sample_param_name = sample_key,
        cell_type = lambda wildcards, output: output[0].split("/")[-2]
    singularity:
        envs['multiome']
    threads:
        64
    resources:
        runtime=240, mem_mb=3000000, disk_mb=500000, slurm_partition='largemem'
    script:
        'scripts/atac_fragment_pseudobulk.py'

rule atac_celltype_call_peaks:
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
        envs['multiome']
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
        envs['scenic']
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        'scripts/MACS_consensus.py'


rule merged_consensus_peak_anndata:
    input:
        consensus_bed = work_dir + '/data/consensus_regions.bed',
        fragment_file=expand(
            data_dir + '{sample}/outs/atac_fragments.tsv.gz',
            sample=samples
            )
    output:
        merged_atac_anndata = work_dir + '/atlas/03_consensus_peak_atac.h5ad',
        output_files = expand(
            data_dir + '{sample}/outs/02_{sample}_anndata_filtered_atac.h5ad',
            sample=samples
            ),
    singularity:
        envs['multiome']
    threads:
        32
    resources:
        runtime=1440, mem_mb=3000000, slurm_partition='largemem'
    script:
        'scripts/atac_merge.py'
    
rule atac_spectral:
    input:
        merged_atac_anndata = work_dir + '/atlas/03_consensus_peak_atac.h5ad'
    output:
        merged_atac_anndata = work_dir + '/atlas/04_modeled_anndata_atac.h5ad'
    params:
        num_features = 100000,
        sample_param = 'sample_id'
    singularity:
        envs['multiome']
    threads:
        32
    resources:
        runtime=1440, mem_mb=250000
    script:
        'scripts/atac_spectral.py'

"""========================================================================="""
"""                               MULTI portion                             """
"""========================================================================="""

rule filter_rna_atac:
    input:
        rna_anndata =data_dir+'{sample}/outs/02_{sample}_anndata_filtered_rna.h5ad',
        atac_anndata = data_dir+'{sample}/outs/02_{sample}_anndata_filtered_atac.h5ad'
    output:
        atac_anndata = data_dir+'{sample}/outs/03_{sample}_anndata_filtered_atac.h5ad',
        rna_anndata = data_dir+'{sample}/outs/03_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['multiome']
    resources:
        runtime=30, mem_mb=50000, slurm_partition='quick'
    script:
        work_dir+'scripts/atac_filter.py'

rule merge_multiome_rna:
    input:
        rna_anndata=expand(
            data_dir+'{sample}/outs/03_{sample}_anndata_filtered_rna.h5ad', 
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'atlas/03_filtered_anndata_rna.h5ad'
    singularity:
        envs['multiome']
    params:
        samples=samples
    resources:
        runtime=120, mem_mb=300000, disk_mb=10000#, slurm_partition='largemem' 
    script:
        work_dir+'scripts/rna_merge.py'

rule merge_multiome_atac:
    input:
        atac_anndata=expand(
            data_dir+'{sample}/outs/03_{sample}_anndata_filtered_atac.h5ad', 
            sample=samples
            )
    output:
        merged_atac_anndata = work_dir+'atlas/03_filtered_anndata_atac.h5ad'
    singularity:
        envs['multiome']
    params:
        samples=samples
    resources:
        runtime=720, mem_mb=3000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_atac.py'

rule merge_multiome:
    input:
        merged_atac_anndata = work_dir+'atlas/04_modeled_anndata_atac.h5ad',
        merged_rna_anndata = work_dir+'atlas/05_QC_filtered_anndata_rna.h5ad'
    output:
        multiome_object = work_dir+'atlas/03_merged_multiome.h5mu'
    singularity:
        envs['multiome']
    params:
        sample_key=sample_key
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        work_dir+'scripts/merge_multiome.py'
        
rule multiome_feature_selection:
    input:
        multiome_object = work_dir+'atlas/03_merged_multiome.h5mu'
    output:
        multiome_object = work_dir+'atlas/04_highly_variable_multiome.h5mu'
    singularity:
        envs['multiome']
    params:
        hvg = 3000,
        hvp = 20000,
        sample_key = sample_key
    resources:
        runtime=480, mem_mb=1500000, slurm_partition='largemem'
    script:
        work_dir+'scripts/multiome_feature_selection.py'

rule multivi:
    input:
        multiome_object = work_dir+'atlas/04_highly_variable_multiome.h5mu'
    output:
        multiome_object = work_dir+'atlas/05_highly_variable_multivi_multiome.h5mu',
        model_history = work_dir+'data/model_elbo/multiome_model_history.csv'
    params:
        model = work_dir+'data/models/multiome_polish/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/multiome_model.sh {input.multiome_object} {params.sample_key} {output.model_history} {output.multiome_object} {params.model}'

rule transfer_UMAP:
    input:
        multiome_object = work_dir+'/atlas/03_merged_multiome.h5mu',
        hvg_multiome_anndata = work_dir + '/atlas/05_highly_variable_multivi_multiome.h5mu'
    output:
        multiome_object = work_dir + '/atlas/06_polished_multiome.h5mu'
    singularity:
        envs['multiome']
    resources:
        runtime=1440, mem_mb=1000000, slurm_partition='largemem'
    script:
        work_dir+'/scripts/multivi_to_UMAP.py'

rule pychromvar:
    input:
        merged_multiome = work_dir + '/atlas/07_annotated_multiome.h5mu',
        reference_genome = reference_genome
    output:
        merged_multiome = work_dir+'atlas/08_multiome.h5mu'
    params:
        chunk_size = 100000
    singularity:
        envs['multiome']
    threads:
        16
    resources:
        runtime=2880, mem_mb=1000000, slurm_partition='largemem'
    script:
        'scripts/pychromvar.py'

"""========================================================================="""
"""                            ANALYSIS portion                             """
"""========================================================================="""

rule rna_pseudobulk:
    input:
        rna_anndata = work_dir + '/atlas/07_polished_anndata_rna.h5ad'
    output:
        pseudo_rna = work_dir+'atlas/pseudobulked_rna.h5ad'
    params:
        sample_key = 'Sample_ID',
        separating_cluster = 'celltype',
        min_cells = 10
    singularity:
        envs['multiome']
    resources:
        runtime=120, mem_mb=250000, slurm_partition='quick'
    script:
        'scripts/rna_pseudobulk.py'

rule atac_pseudobulk:
    input:
        merged_atac_anndata = work_dir + '/atlas/04_modeled_anndata_atac.h5ad'
    output:
        pseudo_atac = work_dir+'atlas/pseudobulked_atac.h5ad'
    params:
        sample_key = 'sample_id',
        separating_cluster = 'celltype',
        min_cells = 10
    singularity:
        envs['multiome']
    resources:
        runtime=120, mem_mb=200000, slurm_partition='quick'
    script:
        'scripts/atac_pseudobulk.py'

rule chromvar_pseudobulk:
    input:
        merged_multiome = work_dir+'atlas/multiome_chromvar_atlas.h5mu'
    output:
        pseudobulk_chromvar = work_dir+'atlas/pseudobulked_chromvar.h5ad'
    params:
        sample_key = 'SampleID',
        separating_cluster = 'celltype'
    singularity:
        envs['multiome']
    resources:
        runtime=120, mem_mb=200000, slurm_partition='quick'
    script:
        'scripts/chromvar_pseudobulk.py'

rule cell_fraction_plot_and_test:
    input:
        merged_rna_anndata = work_dir+'atlas/07_polished_anndata_rna.h5ad'
    output:
        fraction_boxplot = work_dir+'figures/cell_count_by_disease_and_celltype_boxplot.svg',
        corrected_ztest_results = work_dir+'data/celltype_fraction_ztest_results.csv'
    params:
        sample_key = sample_key,
        disease_param = disease_param,
        separating_cluster = 'celltype',
        control = 'control',
        diseases = ['PD', 'LBD'],
        separating_value_dict = dict(zip(['control', 'PD', 'LBD'], ['#7f7f7f', '#5ab4e5', '#d36027']))
    singularity:
        envs['multiome']
    script:
        work_dir + '/scripts/cell_fraction_test_plot.py'

rule cell_cell_communication:
    input:
        merged_rna_anndata = work_dir+'atlas/07_polished_anndata_rna.h5ad',
    output:
        cell_cell_communication_data = work_dir+'data/CCC/combined/CCC_celltype_results.csv'
    params:
        control = control,
        disease_param = disease_param
    threads:
        64
    resources:
        disk_mb=200000, mem_mb=200000, slurm_partition='quick'
    script:
        'scripts/cell_cell_communication.py'

rule DEG:
    input:
        pseudo_rna = work_dir + '/atlas/pseudobulked_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/DGE_Dreampy_results.csv'
    params:
        celltype_params = 'celltype',
        celltypes = cell_types,
        diagnosis_param = disease_param,
        control = control,
        diagnosis_control = [control] + diseases,
        sample_key=sample_key,
        formula = "~ Primary Diagnosis + Age + Sex + (1|Use_batch) + (1|Brain_bank)"
    singularity:
        envs['multiome']
    threads:
        64
    resources:
        runtime=180, mem_mb=200000, slurm_partition='quick'
    script:
        'scripts/rna_DEG.py'

rule differential_cell_cell_communication:
    input:
        merged_rna_anndata = work_dir+'atlas/07_polished_anndata_rna.h5ad',
        merged_DGE_data = work_dir + '/data/significant_genes/rna_unfiltered_gene_hits.csv'
    output:
        differential_cell_cell_communication_data = work_dir + '/data/CCC/differential_CCC_by_{sep_param}_{disease}_pairs.csv'
    params:
        disease = lambda wildcards, output: output[0].split("_")[-2],
        disease_param = disease_param,
        sep_param = lambda wildcards, output: output[0].split("_")[-3],
    script:
        'scripts/rna_differential_cell_cell_communication.py'

rule disease_gsea:
    input:
        adata_path = work_dir+'atlas/07_polished_anndata_rna.h5ad',
        ontologies = work_dir+'input/ontologies.csv'
    output:
        cell_disease_GSEA =  work_dir+'data/GSEA/{separating_cluster}/GSEA_{separating_cluster}_{cell_type}_{control}_{disease}_results.csv'
    params:
        disease_param = disease_param,
        control = lambda wildcards, output: output[0].split("_")[-3],
        separating_cluster = lambda wildcards, output: output[0].split("_")[-5],
        cell_type = lambda wildcards, output: output[0].split("_")[-4],
        disease = lambda wildcards, output: output[0].split("_")[-2]
    singularity:
        envs['multiome']
    threads:
        64
    resources:
        runtime=960, mem_mb=1000000, slurm_partition='largemem' 
    script:
        'scripts/rna_GSEA.py'

rule disease_great:
    input:
        DAR_path =  work_dir+'data/significant_genes/atac/atac_{cell_type}_{control}_{disease}_DAR.csv',
        tss_file =  work_dir+'input/tss_from_great.bed',
        chr_sizes_file =  work_dir+'input/chr_size.bed',
        annotation_file =  work_dir+'input/ontologies.csv',
    output:
        cell_disease_peaks = work_dir+'data/celltypes/{cell_type}/{cell_type}_{control}_{disease}_DAR_peaks.bed',
        cell_disease_GREAT = work_dir+'data/celltypes/{cell_type}/{cell_type}_{control}_{disease}_GREAT_peaks.csv'
    singularity:
        envs['multiome']
    resources:
        runtime=2880
    script:
        'scripts/atac_GREAT.py'

rule celltype_great:
    input:
        consensus_bed = work_dir + '/data/consensus_regions.bed',
        tss_file =  work_dir+'input/tss_from_great.bed',
        chr_sizes_file =  work_dir+'input/chr_size.bed',
        annotation_file =  work_dir+'input/ontologies.csv',
    output:
        cell_disease_GREAT = work_dir+'data/celltypes/{cell_type}/{cell_type}_GREAT_peaks.csv'
    params:
        cell_types = cell_types,
        celltype = lambda wildcards: wildcards.cell_type
    singularity:
        envs['multiome']
    resources:
        runtime=2880
    script:
        'scripts/atac_celltype_GREAT.py'

rule celltype_overlapping_peaks:
    input:
        peak_files = expand(
            work_dir+'data/celltypes/{celltype}/{celltype}_{condition}_peaks.bed',
            condition = diseases + [control],
            allow_missing = True
        )
    output:
        celltype_overlapping_celltype_peaks = work_dir+'data/celltypes/{celltype}/{celltype}_overlapping_peaks.csv'
    singularity:
        envs['atac_fragment']
    resources:
        slurm_partition='quick'
    script:
        'scripts/overlapping_peaks.py'

rule atac_merged_coaccessibilty:
    input:
        celltype_atac = work_dir + '/atlas/multiome_wnn.h5mu'
    output:
        celltype_atac = work_dir + '/atlas/04_coaccessible_anndata_atac.h5ad',
        circe_network = work_dir+'data/circe_network.csv'
    singularity:
        envs['multiome']
    threads:
        16
    resources:
        runtime=1440, mem_mb=1500000, slurm_partition='largemem'
    script:
        'scripts/circe_by_celltype.py'

rule gene_motif_linkage:
    input:
        pseudobulk_chromvar = work_dir+'atlas/pseudobulked_chromvar.h5ad',
        pseudo_rna = work_dir+'atlas/pseudobulked_rna.h5ad'
    output:
        gene_motif_links = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{diagnosis}_gene_motif_link.csv'
    params:
        diagnosis_param = disease_param,
        celltype_param = 'celltype',
        celltype = lambda wildcards: wildcards.cell_type,
        diagnosis = lambda wildcards: wildcards.diagnosis
    singularity:
        envs['multiome']
    resources:
        runtime=180, mem_mb = 50000, slurm_partition = 'quick'
    script:
        'scripts/gene_motif_linkage.py'

rule gene_peak_linkage:
    input:
        pseudobulked_rna = work_dir+'atlas/pseudobulked_rna.h5ad',
        gene_info = gene_info,
        gene_tss = gene_tss,
        atac_files = expand(
            work_dir+'data/celltypes/{celltype}/{celltype}_{condition}_atac.h5ad',
            condition = [control] + diseases,
            allow_missing = True
        ),
        circe_files = expand(
            work_dir+'data/celltypes/{celltype}/{celltype}_{condition}_circe_network.csv',
            condition = [control] + diseases,
            allow_missing = True
        ),
        bed_files = expand(
            work_dir + '/data/celltypes/{celltype}/{celltype}_{condition}_peaks.bed',
            condition = [control] + diseases,
            allow_missing = True
        )
    output:
        gene_peak_linkage = work_dir+'data/celltypes/{celltype}/{celltype}_promoter_coaccessibility.csv'
    params:
        conditions = [control] + diseases,
        celltype = lambda wildcards: wildcards.celltype
    singularity:
        envs['multiome']
    resources:
        slurm_partition='quick'
    script:
        'scripts/gene_peak_linkage.py'

"""========================================================================="""
"""                            CELLTYPE portion                             """
"""========================================================================="""

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
        merged_rna_anndata = work_dir+'atlas/07_polished_anndata_rna.h5ad',
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_annotated_peaks.bed',
        fragment_files=expand(
            data_dir+'{sample}/outs/atac_fragments.tsv.gz',
            sample=samples,
            )
    output:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/atac.h5ad'
    singularity:
        envs['multiome']
    params:
        pseudobulk_param = 'celltype',
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

rule atac_coaccessibilty:
    input:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/atac.h5ad'
    output:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/atac_circe.h5ad',
        circe_network = work_dir+'data/celltypes/{cell_type}/circe_network_{cell_type}.csv'
    params:
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    singularity:
        envs['multiome']
    threads:
        16
    resources:
        runtime=1440, mem_mb=1500000, slurm_partition='largemem'
    script:
        'scripts/circe_by_celltype.py'

rule fragments_pseudobulk_cell_disease:
    input:
        merged_rna_anndata = work_dir+'atlas/07_polished_anndata_rna.h5ad',
        fragment_file=expand(
            data_dir+'{sample}/outs/atac_fragments.tsv.gz',
            sample=samples,
            )
    output:
        pseudo_fragment_files = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_fragments.bed'
    params:
        pseudobulk_param = 'cell_type',
        samples=samples,
        sample_param_name = sample_key,
        cell_type = lambda wildcards: wildcards.cell_type,
        disease = lambda wildcards: wildcards.disease,
        disease_param = disease_param
    singularity:
        envs['multiome']
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
        envs['multiome']
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
        blacklist = work_dir + '/input/hg38-blacklist.bed'
    singularity:
        envs['multiome']
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
        merged_rna_anndata = work_dir+'atlas/07_polished_anndata_rna.h5ad',
        cell_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_peaks.bed',
        cell_annotated_bedfile = work_dir + '/data/celltypes/{cell_type}/{cell_type}_{disease}_annotated_peaks.bed',
        fragment_files=expand(
            data_dir+'{sample}/outs/atac_fragments.tsv.gz',
            sample=samples,
            )
    output:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_atac.h5ad'
    singularity:
        envs['multiome']
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
        celltype_atac = work_dir+'data/celltypes/{cell_type}/{cell_type}_{control}_{disease}_atac.h5ad'
    output:
        celltype_atac = work_dir+'data/celltypes/{cell_type}/{cell_type}_{control}_{disease}_atac_circe.h5ad',
        circe_network = work_dir+'data/celltypes/{cell_type}/{cell_type}_{control}_{disease}_circe_network.csv'
    params:
        cell_type = lambda wildcards: wildcards.cell_type
    singularity:
        envs['multiome']
    threads:
        8
    resources:
        runtime=600, mem_mb=400000, slurm_partition='largemem'
    script:
        'scripts/circe_by_celltype.py'

"""========================================================================="""
"""                              MOTIF portion                              """
"""========================================================================="""

rule motif_enrichment:
    input:
        atac_anndata = work_dir+'atlas/04_modeled_anndata_atac.h5ad',
        ref_genome = reference_genome,
        TF_motifs = motifs
    output:
        motif_enrichment = work_dir+'data/motif_enrichment.csv'
    params:
        control = control,
        cell_type = 'celltype',
        disease_param = disease_param
    singularity:
        envs['multiome']
    resources:
        runtime=240, disk_mb=300000, mem_mb=200000
    script:
        'scripts/atac_motif_enrichment.py'

rule differential_motif_enrichment:
    input:
        output_DAR_data = work_dir+'data/DARs/{separating_cluster}/DAR_{separating_cluster}_{cell_type}_{control}_{disease}_DAR.csv',
        cell_type_atac = work_dir+'data/celltypes/{cell_type}/atac.h5ad',
        TF_motifs = motifs,
        ref_genome = reference_genome
    output:
        differential_motif_dataframe = work_dir+'data/DMEs/{separating_cluster}/DME_{separating_cluster}_{cell_type}_{control}_{disease}_results.csv'
    params:
        disease_param = disease_param,
        design_factors = design_covariates,
        control = lambda wildcards, output: output[0].split("_")[-3],
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-4],
        separating_cluster = lambda wildcards, output: output[0].split("_")[-5],
    singularity:
        envs['multiome']
    resources:
        runtime=240, disk_mb=300000, mem_mb=200000
    script:
        'scripts/differential_motif_enrichment.py'

rule barcode_merge:
    input:
        cell_annotate = work_dir+'data/rna_cell_annot.csv',
        metadata_table = metadata_table
    output:
        annotate_metadata_table = work_dir+'data/barcode_cell_annotation.csv'
    params:
        disease_param = disease_param,
        sample_key = sample_key
    resources:
        slurm_partition='quick' 
    run:
        import pandas as pd
        metadata_df = pd.read_csv(input.metadata_table)
        sample_disease = dict(zip(metadata_df[params.sample_key], metadata_df[params.disease_param]))
        cell_barcodes = pd.read_csv(input.cell_annotate)
        cell_barcodes['sample'] = ['_'.join(x.split('_')[1:]) for x in cell_barcodes['atlas_identifier']]
        cell_barcodes['disease'] = [sample_disease[x] for x in cell_barcodes['sample']]
        cell_barcodes['barcode'] = [x.split('_')[0] for x in cell_barcodes['atlas_identifier']]
        cell_barcodes.to_csv(output.annotate_metadata_table, index=False)

rule filter_celltype_condition_samples:
    input:
        annotate_metadata_table = work_dir+'data/barcode_cell_annotation.csv'
    output:
        batch_sample_celltype_disease_df = work_dir+'data/batch_sample_celltype_disease.csv'
    params:
        seq_batch_key = seq_batch_key
    run:
        import pandas as pd
        annotate_metadata_table = pd.read_csv(input.annotate_metadata_table)
        combindation_df = annotate_metadata_table.groupby([params.seq_batch_key, 'sample', 'cell_type', 'disease']).count().reset_index()
        bscd_df = combindation_df[combindation_df['barcode'] != 0][[params.seq_batch_key, 'sample', 'cell_type', 'disease']]
        bscd_df.to_csv(output.batch_sample_celltype_disease_df, index=False)

rule barcode_filter:
    input:
        annotate_metadata_table = work_dir+'data/barcode_cell_annotation.csv'
    output:
        cell_disease_barcodes = temp(work_dir+'data/celltypes/{cell_type}/batch{batch}_{sample}_{cell_type}_{disease}_barcodes.txt')
    resources:
        slurm_partition='quick'
    shell:
        "python scripts/filter_barcode.py {input.annotate_metadata_table} {wildcards.cell_type} {wildcards.sample} {output.cell_disease_barcodes}"

def filter_celltype_condition_samples_seq(wildcards):
    df = pd.read_csv(work_dir+'data/barcode_cell_annotation.csv')
    df = df[(df['celltype'] == wildcards.cell_type) & (df['disease'] == wildcards.disease)][['celltype', seq_batch_key, 'sample', 'disease']].drop_duplicates()
    return [data_dir + f"batch{str(df.loc[x, 'Use_batch'])}/Multiome/{str(df.loc[x, 'sample'])}-ARC/outs/atac_{str(df.loc[x, 'celltype'])}_{str(df.loc[x, 'disease'])}.bam" for x in df.index]

rule celltype_sample_filter_bam:
    input:
        cell_disease_barcodes = work_dir+'data/celltypes/{cell_type}/batch{batch}_{sample}_{cell_type}_{disease}_barcodes.txt',
        input_bam = data_dir+'{sample}/outs/atac_possorted_bam.bam'
    output:
        sample_filter_bam = data_dir+"{sample}outs/atac_{cell_type}_{disease}.bam",
        output_header = temp(data_dir+"{sample}/outs/atac_{cell_type}_{disease}_header"),
        output_body = temp(data_dir+"{sample}/outs/atac_{cell_type}_{disease}_body.sam"),
        output_sam = temp(data_dir+"{sample}/outs/atac_{cell_type}_{disease}.sam")
    singularity:
        envs['multiome']
    threads:
        16
    resources:
        slurm_partition='quick'
    shell:
        "samtools view -H {input.input_bam} -@ 16 > {output.output_header} \n"
        "samtools view {input.input_bam} -@ 16 | LC_ALL=C grep -F -f {input.cell_disease_barcodes} > {output.output_body} \n"
        "cat {output.output_header} {output.output_body} > {output.output_sam} \n"
        "samtools view -b {output.output_sam} -@ 16 > {output.sample_filter_bam}"

rule pseudobulk_bams:
    input:
        sample_filter_bam = filter_celltype_condition_samples_seq
    output:
        sorted_pseudobulk_bam = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}.bam',
        presorted_pseudobulk_bam = temp(work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_presort.bam')
    params:
        batch_sample_celltype_disease_df = work_dir+'data/batch_sample_celltype_disease.csv'
    singularity:
        envs['multiome']
    threads:
        16
    resources:
        runtime=960, mem_mb=300000
    shell:
        "samtools cat {input.sample_filter_bam} -o {output.presorted_pseudobulk_bam} \n"
        "samtools sort {output.presorted_pseudobulk_bam} -@ 16 -o {output.sorted_pseudobulk_bam}"

rule celltype_disease_ATACorrect:
    input:
        sorted_pseudobulk_bam = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}.bam',
        blacklist = work_dir + '/input/hg38-blacklist.bed',
        cell_type_peaks = work_dir+'data/celltypes/{cell_type}/{cell_type}_peaks.bed',
        ref_genome = reference_genome
    output:
        corrected_bigwig = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_corrected.bw'
    params:
        ATACorrect_outdir = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect',
        prefix = '{cell_type}_{disease}'
    singularity:
        envs['tobias']
    threads:
        16
    resources:
        runtime=960, mem_mb=300000
    shell:
        'TOBIAS ATACorrect --bam {input.sorted_pseudobulk_bam} --genome {input.ref_genome} --blacklist {input.blacklist} --peaks {input.cell_type_peaks} --outdir {params.ATACorrect_outdir} --prefix {params.prefix} --cores {threads}'

rule celltype_disease_score_bigwig:
    input:
        corrected_bigwig = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_corrected.bw',
        regions = work_dir+'data/celltypes/{cell_type}/{cell_type}_peaks.bed',
    output:
        footprinted_bigwig = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_footprints.bw'
    singularity:
        envs['tobias']
    threads:
        64
    resources:
        runtime=960, mem_mb=300000
    shell:
        'TOBIAS FootprintScores --signal {input.corrected_bigwig} --regions {input.regions} --output {output.footprinted_bigwig} --cores {threads}'

rule control_comparison_score_bigwig:
    input:
        corrected_bigwig = work_dir+'data/celltypes/{cell_type}/{cell_type}_control_ATACorrect/{cell_type}_control_corrected.bw',
        regions = work_dir+'data/consensus_regions.bed'
    output:
        control_footprint_bigwig = work_dir+'data/celltypes/{cell_type}/{cell_type}_control_ATACorrect/{cell_type}_control_comparison_footprints.bw'
    singularity:
        envs['tobias']
    threads:
        64
    resources:
        runtime=960, mem_mb=300000
    shell:
        'TOBIAS FootprintScores --signal {input.corrected_bigwig} --regions {input.regions} --output {output.control_footprint_bigwig} --cores {threads}'

rule disease_footprinting:
    input:
        motifs = motifs,
        control_bw = work_dir+'data/celltypes/{cell_type}/{cell_type}_{control}_ATACorrect/{cell_type}_{control}_corrected.bw',
        disease_bw = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_ATACorrect/{cell_type}_{disease}_corrected.bw',
        peaks=work_dir+'data/celltypes/{celltype}/{celltype}_overlapping_peaks.bed',
        genome=reference_genome
    output:
        control_disease_motif_data = work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_{control}_BINDetect/bindetect_results.txt'
    params:
        outdir=work_dir+'data/celltypes/{cell_type}/{cell_type}_{disease}_{control}_BINDetect'
    singularity:
        envs['tobias']
    threads:
        16
    resources:
        runtime=180, mem_mb=200000, slurm_partition='quick'
    shell:
        'TOBIAS BINDetect --motifs {input.motifs} --signals {input.control_bw} {input.disease_bw} --genome {input.genome} --peaks {input.peaks}  --outdir {params.outdir} --cores {threads}'

rule control_footprinting:
    input:
        motifs = motifs,
        control_bw = work_dir+'data/celltypes/{cell_type}/{cell_type}_control_ATACorrect/{cell_type}_control_comparison_footprints.bw',
        peaks=work_dir+'data/celltypes/{celltype}/{celltype}_overlapping_peaks.bed',
        genome=reference_genome
    output:
        control_motif_data = work_dir+'data/celltypes/{cell_type}/{cell_type}_control_BINDetect/bindetect_results.txt'
    params:
        outdir=work_dir+'data/celltypes/{cell_type}/{cell_type}_control_BINDetect'
    singularity:
        envs['tobias']
    threads:
        16
    resources:
        runtime=180, mem_mb=200000, slurm_partition='quick'
    shell:
        'TOBIAS BINDetect --motifs {input.motifs} --signals {input.control_bw}  --genome {input.genome} --peaks {input.peaks}  --outdir {params.outdir} --cores {threads}'
