## scMAVERICS - snRNAseq processing

This tutorial runs through taking the output of Cellbender (which should be run outside of a notebook) for a handful of samples and integrate the samples into a cell type annotated object. If the dataset in question is only transcriptomics, not multiomics, this analysis should be sufficent.

This section borrows heavily from the [scanpy](https://scanpy.scverse.org/en/1.12.x/) (version 1.11.5) introduction tutorial, with modifications made for incorporating packages
- [scVI](https://scvi-tools.org/) (version 1.4.2)
- [decoupler](https://decoupler.scverse.org/en/stable/) (version 2.1.6)
- [scrublet](https://github.com/swolock/scrublet)

## Package and data loading


```python
import pandas as pd
import scanpy as sc
import decoupler as dc
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
```

### Define working and data directories


```python
data_dir = '/data/CARD_singlecell/Brain_atlas/SN_multiome/'
work_dir = '/data/CARD_singlecell/SN_atlas/' 
```

### Load metadata

As stated before, metadata should be 


```python
metadata_df = pd.read_csv(work_dir + 'input/metadata.csv')
metadata_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Sample_ID</th>
      <th>Use_batch</th>
      <th>Age</th>
      <th>PMI</th>
      <th>Sex</th>
      <th>Brain_bank</th>
      <th>Primary Diagnosis</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>6178</td>
      <td>6</td>
      <td>68</td>
      <td>22.5</td>
      <td>Male</td>
      <td>Human Brain and Spinal Fluid Resource Center</td>
      <td>control</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3284</td>
      <td>1</td>
      <td>78</td>
      <td>11.5</td>
      <td>Female</td>
      <td>Human Brain and Spinal Fluid Resource Center</td>
      <td>control</td>
    </tr>
    <tr>
      <th>2</th>
      <td>6384</td>
      <td>3</td>
      <td>50</td>
      <td>22.3</td>
      <td>Female</td>
      <td>Human Brain and Spinal Fluid Resource Center</td>
      <td>PD</td>
    </tr>
    <tr>
      <th>3</th>
      <td>5638</td>
      <td>3</td>
      <td>77</td>
      <td>19.3</td>
      <td>Male</td>
      <td>Human Brain and Spinal Fluid Resource Center</td>
      <td>PD</td>
    </tr>
  </tbody>
</table>
</div>



### Define parameters


```python
# Define keys
sample_key = 'Sample_ID'
batch_key = 'Use_batch'

# Define list of samples
samples = metadata_df[sample_key]
batches = metadata_df[batch_key]
```

# Sample preprocessing

This analysis begins after running cellbender. As this process is GPU/CPU intensive and takes a long time to run, it is best done via terminal with the command:
```
cd path/to/sample/outs/

cellbender remove-background \
--input raw_feature_bc_matrix.h5 \
--output cellbender_gex_counts_filtered.h5 \
--checkpoint ckpt.tar.gz --cuda --posterior-batch-size 64
```

This will run for a while for each sample. Make sure you are in the desired sample folder when running, as the `ckpt.tar.gz` file will be created in the current working directory, and checkpoint files will overwrite one another, leading to errors and stalled Cellbender runs.

## Preprocessing



### rna_preprocess




```python
# Define input files
step_1_files = [data_dir + f'/batch{batches[i]}/Multiome/{samples[i]}-ARC/outs/cellbender_gex_counts_filtered.h5' for i in range(len(samples))]
# Define the output files (corresponding to sample)
step_2_adatas = {}
```


```python

```

### rna_filter




```python

```

### rna_merge




```python

```

### rna_feature_selection




```python

```

### rna_model




```python

```

### rna_latent_transfer




```python

```

### rna_annotate




```python

```

### rna_cluster_based_QC




```python

```
