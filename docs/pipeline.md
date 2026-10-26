# scMAVERICS

**Using the snakemake pipeline**

# Installation and setup

## Downloading the pipeline

### Snakemake pipeline
Begin by cloning the scMAVERICS repo with
```
git clone https://github.com/NIH-CARD/scMAVERICS.git
cd scMAVERICS
```

There is no need to worry about python environments. All scripts in the Snakemake-enabled pipeline use singularity ([apptainer](https://apptainer.org/documentation/)) and [snakemake](https://snakemake.readthedocs.io/en/v8.0.0/), with singularity images (`.sif` files) being pulled from [Quay.io](https://quay.io/). Each step of the analysis assumes this is being run on a Slurm-enabled HPC, where GPUs and large-memory nodes are available. 

Large memory nodes are required for handling large (>50) sample multiome datasets and GPUs are required for running samples in a timely manner. 

### Local environment
However, if you don't have access to an HPC or have not enough samples to justify running the entire pipeline, the environment used for most steps can be recreated with:

```
git clone https://github.com/NIH-CARD/scMAVERICS.git
cd scMAVERICS
conda create -n multiome
conda activate multiome
pip install -r requirements.txt
```

This environment is used in example analysis in Jupyter notebooks, beginning with [snRNAseq analysis](https://scmaverics.readthedocs.io/en/latest/snRNAseq/). 


## Importing data

### Starting data

Sample data should be stored in separate folders corresponding to sample. This analysis assumes multiome runs of [10X-based single-nuclei multiome sequencing](https://www.10xgenomics.com/support/epi-multiome) has been performed. Subsequently, data processing should be done before running this pipeline, ideally with [CellRanger-ARC](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest). 

Additionally, it is assumed that the output folders of CellRanger-ARC correspond to a `sample_key` column in the metadata, and the folders contain:

```
cellranger_arc_outputs
├── Sample_1
│   └──outs
│       ├── raw_feature_bc_matrix.h5
│       ├── atac_possorted_bam.bam
│       ├── atac_fragments.tsv.gz
│       └── atac_fragments.tsv.gz.tbi
├── Sample_2
│   └──outs
│       ├── raw_feature_bc_matrix.h5
│       ├── atac_possorted_bam.bam
│       ├── atac_fragments.tsv.gz
│       └── atac_fragments.tsv.gz.tbi
...
```

### Metadata and reference files



## Config file setup

Config files ensure consitency and reproducibility. In this case, 

# snRNAseq steps

### Cellbender

### Preprocessing and integration

### scVI and cell type annotation

# snATACseq steps

### Preprocessing and integration

### Cell type and consensus peaks

## Multiome object creation

### Merging modalities and MultiVI

### TF motif inference with pyChromvar

# Differential analysis

### Cell-cell communication

### Differentially expressed genes

### Differentially enriched motifs

### Differential cell-cell communication

### Gene/peak-set enrichment

### Peak co-accessibility

### Gene-motif/peak linkage

# ATAC fragment analysis

### Cell type and diagnosis peaks/bigwigs

### Cell type and diagnosis peak co-accessibility
