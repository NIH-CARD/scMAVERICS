# scMAVERICS

TL;DR 

Clone the repo, add metadata, update config.yaml, select the step you want completed, run `bash snakemake.sh`.

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

Snakemake doesn't know where your data lives! Preprocessing and aggregating individual samples is handled via generated paths, dependent on a *metadata* file provided by the user.

**🚨 CRITICAL 🚨**

No metadata is provided, it is critical the user provides appropriate metadata as a file `metadata.csv` in the `input` folder. Any metadata entered in file will be added to the final object, but the bare minimum there should be columns of `Sample_ID`, `diagnosis`, and `sex`.  These parameters will be used in downstream analysis.

```
file: input/metadata.csv ⬇

| Sample_ID   | diagnosis   | sex            | add. metadata
|-------------|-------------|----------------|-----
| Sample_1    | Control     | male           | ...
| Sample_2    | AD          | female         | ...
```
## Config file setup

Config files ensure consitency and reproducibility. In this case, the config file is what truly differentiates one analysis from another. This includes:
```
- Data and working directories
- Desired output files
- System parameters
- Reference files (e.g. reference genomes)
- Parameters of the dataset
- QC values
```

**🚨 CRITICAL 🚨**

The config file is unique to the input dataset. Only the scientists working with the data know which reference file was for Cellranger, or what QC thresholds make sense for the data. Sometimes parameters need to be changed post-hoc, in which case it is import that this file is version controlled somewhere, either by a fork from the original repo or a backup copy.

### Data and working directories
Both `data_dir` and `work_dir` need to be updated to reflect where the outputs of CellRanger-ARC are stored (`data_dir`) and where what is the location of this pipeline (`work_dir`).

### Output files
Snakemake works by defined the desired output file and determining what rules need to be executed in order to yield the file. From that, the file specified in `final_output` simplifies determining what to call for the output. A list of possible outputs is outlined above `final_output`. It is possible enter the final object `'atlas/08_multiome.h5mu'`, but it is recommended that each step be run sequentially, in order to inspect the intermediary datasets.

### System parameters
It is assumed that this pipeline is being run on a HPC cluster powered by Slurm, but the limits of individual partitions vary. Specify the parameters of the user's node options here, between `short_gpu`, `multi_gpu`, `quick`, `long`, `large_1`, and `large_2`.

### Reference files
It is assumed that when CellRanger-ARC was run, there were reference files used to map reads. `reference_genome`, `genome_length`, `gene_tss`, and `gene_info` should already be downloaded from running CellRanger-ARC, simply specify their paths here.

Additionally, other input files are required. `metadata`, as mentioned above, is critical for running this pipeline. Additionally, `gene_list` is required for annotating cell types and `motif_file` is required for TFBS activity inference (download in .meme format from [here](https://jaspar.elixir.no/)).

### Sample parameters
The parameters `sample_key`, `sex_key`, and `disease_param` should be specified as columns within the metadata file. Additionally, `control` should be the name of the reference sample state and `diagnoses` should be a list of altered sample states. `covariates` should include the fixed parameters that co-vary between donor states.

### QC thresholds
Values here should be determined dataset-to-dataset, with suggested values provided. 