# scMAVERICS


**Using the snakemake pipeline**

## Installation

### Snakemake pipeline
Begin by cloning the scMAVERICS repo with
```
git clone https://github.com/NIH-CARD/scMAVERICS.git
cd scMAVERICS
```

There is no need to worry about python environments. All scripts in the Snakemake-enabled pipeline use singularity (apptainer) and snakemake, with singularity images (`.sif` files) being pulled from Quay.io. Each step of the analysis assumes this is being run on a Slurm-enabled HPC, where GPUs and large-memory nodes are available. 

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

This environment is used in the example notebook pages.

## Importing data

## Config file setup

## snRNAseq steps

### Cellbender

### Preprocessing and integration

### scVI and cell type annotation

## snATACseq steps

### Preprocessing and integration

### Cell type and consensus peaks

## Multiome object creation

### Merging modalities and MultiVI

### TF motif inference with pyChromvar

## Differential analysis

### Cell-cell communication

### Differentially expressed genes

### Differentially enriched motifs

### Differential cell-cell communication

### Gene/peak-set enrichment

### Peak co-accessibility

### Gene-motif/peak linkage

## ATAC fragment analysis

### Cell type and diagnosis peaks/bigwigs

### Cell type and diagnosis peak co-accessibility
