# scMAVERICS


## single-cell Multiome Analysis using Variational-inference and Enhancer-driven Regulatory-networks to Inform Cell-atlas Structure

Focus on the analysis, not the processing. 

The workflow was built to minimize the time required to process data and get the single-cell scientist analyzing high-quality processed data and generating figures. We tried to strike a balance between giving the user control. The scripts and Snakemake pipeline work for the Multiome datasets output from the 10X Multiome kit (https://www.10xgenomics.com/support/epi-multiome), generated by the Single Cell Expert group at Center for Alzheimer's and Related Dementias at the NIH. The multiome pipeline was built to take several batches of human brain single-nuclei sequencing samples and process them into a multiome atlas object for further analysis. 
 
The modules of Scanpy (https://github.com/scverse/scanpy), SCVI (https://github.com/scverse/scvi-tools), snapATAC2 (https://github.com/kaizhang/SnapATAC2), pycisTopic (https://github.com/aertslab/pycisTopic), and MACS (https://github.com/macs3-project/MACS) are utilized heavily to produce a multiome atlas with minimal batch effects. 

## Pipeline

![screenshot](images/scMAVERICS.png)

### To get started

Copy this repository to where you will be working with your data. This folder will be where output data is stored, while intermediary files will be stored in a separate folder to be defined by the user. It is important that the output of your CellRanger-ARC run has the format: `<data_dir>/<sample>/`!

#### Required inputs:
- Metadata file in .csv format, example in `input/example_metadata.csv`. A minimal metadata file should include:
  - Sequencing batch, called Use_batch in example (indicating the separate sequencing folder for each run of `CellRanger-ARC`, denoted batch"Sequencing batch>")
  - Sample ID, called Sample in the example, indicates where name of each sample
  - Sample comparison, usually a disease or diagnosis, called 'Primary Diagnosis' in the example, as well for the values of this column;
    - Names of the control condition, called 'control' in the example
    - A list of conditions to compare to the control, called 'diseases' in the example
- Cell-typing table with marker genes, in .csv format, see example in `input/example_marker_genes.csv`
- Directory of the CellRanger output, the data here is the output of `CellRanger`, where each sample is contained in a separate data directory folder; starting with `<data_dir>/<sample>/raw_feature_bc_matrix.h5` for each sample

In addition, the `snakefile` requires modifications to fit your project. The top section "Parameter" should be modified for your dataset, include quality control values, where the input metadata and cell/cell gene marker files are stored. Input files should have their values match the parameters section.

#### Outputs:
- RNA and ATAC Multiome atlas object (multiome_atlas.h5ad)
- List of differentially expressed genes and accessible regions (data/significant_genes/(rna or atac)/(celltype)_(disease)_(DGE or DAR).csv
- Preprocessed and QC-filtered AnnData objects for each sample

#### Current version:
- Uses Singularity images for reproducible runs (scVI modeling will be slow until a GPU-enabled image is created)
- Snakemake runs steps until all output files are created
- Genes used for celltyping are input from `input/example_marker_genes.csv`
- Both RNA and ATAC processing has to be done to run this pipeline 
- Data needs to be stored in a specific heirarchy
- Cellbender needs to be run after CellRanger to compensate for ambient RNA
- Differential Gene Expression and Differential Accessibility of Regions analysis are done 

Once set up, this complete pipeline can be run by simply typing `bash snakemake.sh` in terminal in an HPC running Slurm. This is a work in progress and has not been tested on other devices. 

### Ambient RNA correction (rule cellbender)

Cellbender is used to correct for ambient RNA in the unfiltered RNA .h5 output of CellRanger. This can be run without GPUs, but it is incredibly slow. Current configuration allows for using Biowulf GPUs.

### RNA processing (rule preprocess) 

Transcriptomic data from CellRanger-ARC-2.0 ('''cellbender_gex_counts_filtered.h5''') is read in and processed with Scanpy. QC metrics of percent mitochondria/ribosomal RNA, doublet probability, and cell cycle.

### RNA QC (rule filter_rna) 

Parameters from the processing step are used to filter the cells from each samples based on percent mitochondrial transcripts, probability of being a doublet, and the minimum number of genes observed per cell.

### Individual RNA sample merging into atlas (rule merge_filtered_rna)

Each individual RNA AnnData object are merged into a single QC-filtered object for downstream analysis. This isn't required to be run in a normal workflow.

### ATAC processing (rule atac_preprocess)

ATAC fragment data is converted into an AnnData object with bins used as the measured variable in each cell. One object is created for each sample.

### ATAC QC (rule filter_atac) 

Cells in each sample's ATAC object are filtered for a minimum number of bins per cell. 

### Filtering RNA and ATAC data (rule filter_atac) 

Each sample's QC-filtered RNA and ATAC AnnData objects are filtered for the same cells observed in both samples. Final AnnData objects are saved with a '''03_''' prefix.

### Individual RNA sample merging into atlas (rule merge_multiome_rna)

ATAC-filtered RNA samples are merged from rule filter_atac prior to modeling.

### RNA modeling (rule rna_model) 

Filtered RNA samples are merged into an atlas and multidimensional scaling is performed. A copy of the atlas is made with mitochondiral and ribosomal transcripts removed and only the most variable genes kept. SCVI is used to model the embed dimensions of the atlas, with batch correction, followed by KNN, leiden clustering, and UMAP scaling.

### Cell-typing (rule annotate) 

Cell types of the modeled and clustered RNA atlas are estimated using over-representation analysis and a currated list of cell gene markers.

### Fragment pseudobulk (rule cistopic_pseudobulk)

Use the cell annotated barcodes from the RNA atlas annotation to add fragments into all fragments for all cell types

### Peak calling by cell type (rule MACS2_peak_call)

Call peaks using MACS2 for each cell type, based on pseudobulked fragments

### Consensus peak calling (rule consensus_peaks)

Create consensus peaks from cell type peaks

### Peak by cell AnnData creation (rule cistopic_create_objects)

Populate AnnData objects with fragments contained by peak, by sample

### Create peak atlas (rule cistopic_merge_objects)

Concatenate peak AnnData objects into one object

### Model ATAC peaks (rule atac_peaks_model)

Use poissonVI to model ATAC peaks

### Merging to one multiome object (rule multiome_output) 

Both atlases are merged into a single muon AnnData object for portability.

### Separate atlas into individual celltypes (rule export_celltypes) 

Slice the multiome atlas into individual RNA and ATAC AnnData objects by celltype

### Differential Gene Expression and Differentially Accessible Chromatin (rule DGE and DAR) 
From the individual RNA and ATAC AnnData objects, separated by celltype, pseudobulk and compare the expression of RNA transcripts and ATAC availability, export lists of genes/genome bins with fold-changes and p-values, as well as pseudobulked objects
