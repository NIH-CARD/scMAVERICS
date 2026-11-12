## scMAVERICS - multiome integration

The word integration is used a lot in single cell analysis. In this context, multiome integration really means the combining of gene expression and chromatin accessibility modalities across all samples in a dataset. The resulting object at the end of this processing step will be the final object from which all multiome analyses can be based. For analysis specific to chromatin accessibility peaksets unique to cell types and/or diagnosis state binning, proceed to the next section [Peak calling](peaks.md).

**🚨 HALT! 🚨**

This section assumes you have completed the previous section, [snRNAseq preprocessing](snRNAseq.md) and [snATACseq preprocessing](snATACseq.md), with the object `06_QC_filtered_anndata_rna.h5ad` and `04_modeled_anndata_atac.h5ad` in the `atlas` folder. If not, go back and complete the analysis, preferably with the Snakemake pipeline.

## Package and data loading


```python
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import snapatac2 as snap
import scipy.sparse as sp
import scvi
import pychromvar as pc
from pyjaspar import jaspardb
import scipy
import decoupler as dc
import circe as ci
import matplotlib.pyplot as plt
```

### Define working and data directories


```python
data_dir = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/'
work_dir = '/data/CARD_singlecell/SN_atlas/' 
```

### Load metadata

As stated before, metadata should be imported, for visualizing parameters, 


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




```python
# Define keys
sample_key = 'Sample_ID'
batch_key = 'Use_batch'

# Define list of samples
samples = metadata_df[sample_key].astype(str)
batches = metadata_df[batch_key]
# Ensure sample names are strings
metadata_df[sample_key] = metadata_df[sample_key].astype(str)

# Latent space parameters
high_var_gene_num = 2000
high_var_peak_num = 20000
random_number_seed = 0
max_epoch = 1000
machine_type = 'cpu' # mps if you are running this on a Mac with a metal cpu, gpu on a GPU node

# Reference genome
ref_genome = '/fdb/cellranger-arc/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa' 
```

## Creating the multiome object

### multiome_merge

MuData objects are similar to AnnData objects, where each modality is stored as a `mod` within the object. 


```python
# Read in RNA dataset
rna = sc.read_h5ad(work_dir + '/atlas/06_QC_filtered_anndata_rna.h5ad')

# Read in ATAC dataset
atac = sc.read_h5ad(work_dir + '/atlas/04_modeled_anndata_atac.h5ad')
```


```python
# Create the anndata object
mdata = mu.MuData({"rna": rna, "atac": atac})
mdata
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1471: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("var", axis=0, join_common=join_common)
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1322: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("obs", axis=1, join_common=join_common)





<pre>MuData object with n_obs × n_vars = 37281 × 645578
  2 modalities
    rna:	23019 × 36601
      obs:	&#x27;n_genes_by_counts&#x27;, &#x27;total_counts&#x27;, &#x27;total_counts_rb&#x27;, &#x27;pct_counts_rb&#x27;, &#x27;total_counts_mt&#x27;, &#x27;pct_counts_mt&#x27;, &#x27;Sample_ID&#x27;, &#x27;Use_batch&#x27;, &#x27;Age&#x27;, &#x27;PMI&#x27;, &#x27;Sex&#x27;, &#x27;Brain_bank&#x27;, &#x27;Primary Diagnosis&#x27;, &#x27;S_score&#x27;, &#x27;G2M_score&#x27;, &#x27;phase&#x27;, &#x27;doublet_score&#x27;, &#x27;predicted_doublet&#x27;, &#x27;leiden&#x27;, &#x27;celltype&#x27;
      var:	&#x27;gene_ids&#x27;, &#x27;mt&#x27;, &#x27;rb&#x27;, &#x27;highly_variable&#x27;, &#x27;highly_variable_rank&#x27;, &#x27;means&#x27;, &#x27;variances&#x27;, &#x27;variances_norm&#x27;
      uns:	&#x27;celltype_colors&#x27;, &#x27;hvg&#x27;, &#x27;leiden&#x27;, &#x27;leiden_colors&#x27;, &#x27;neighbors&#x27;, &#x27;umap&#x27;
      obsm:	&#x27;X_scvi&#x27;, &#x27;X_umap&#x27;, &#x27;padj_ulm&#x27;, &#x27;score_ulm&#x27;
      layers:	&#x27;counts&#x27;, &#x27;cpm&#x27;, &#x27;log-norm&#x27;
      obsp:	&#x27;connectivities&#x27;, &#x27;distances&#x27;
    atac:	34934 × 608977
      obs:	&#x27;sample&#x27;, &#x27;n_fragment&#x27;, &#x27;tsse&#x27;
      var:	&#x27;count&#x27;, &#x27;selected&#x27;
      uns:	&#x27;spectral_eigenvalue&#x27;
      obsm:	&#x27;X_spectral&#x27;, &#x27;X_umap&#x27;
      obsp:	&#x27;distances&#x27;</pre>



We can see here that as long as the indices of both modalities `rna` and `atac` have the same convention, the `MuData` object will contain both under `mdata.mod['rna']` and `mdata.mod['atac']`. 

In the previous sections we set up the samples with data 


```python
# Create sample key for batch-corrected MultiVI
mdata.obs[sample_key] = [(x.split('-1_')[1]) for x in mdata.obs_names]
mdata.obs
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
      <th>rna:n_genes_by_counts</th>
      <th>rna:total_counts</th>
      <th>rna:total_counts_rb</th>
      <th>rna:pct_counts_rb</th>
      <th>rna:total_counts_mt</th>
      <th>rna:pct_counts_mt</th>
      <th>rna:Sample_ID</th>
      <th>rna:Use_batch</th>
      <th>rna:Age</th>
      <th>rna:PMI</th>
      <th>...</th>
      <th>rna:G2M_score</th>
      <th>rna:phase</th>
      <th>rna:doublet_score</th>
      <th>rna:predicted_doublet</th>
      <th>rna:leiden</th>
      <th>rna:celltype</th>
      <th>atac:sample</th>
      <th>atac:n_fragment</th>
      <th>atac:tsse</th>
      <th>Sample_ID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CTCCAACTCTTGCAAA-1_6178</th>
      <td>11969</td>
      <td>107576</td>
      <td>3301</td>
      <td>3.068528</td>
      <td>5262</td>
      <td>4.891426</td>
      <td>6178</td>
      <td>6</td>
      <td>68</td>
      <td>22.5</td>
      <td>...</td>
      <td>-1.098352</td>
      <td>G1</td>
      <td>0.089721</td>
      <td>False</td>
      <td>0</td>
      <td>EC</td>
      <td>6178</td>
      <td>343949</td>
      <td>17.177150</td>
      <td>6178</td>
    </tr>
    <tr>
      <th>ACAGGATCACAACCTA-1_6178</th>
      <td>7406</td>
      <td>25161</td>
      <td>222</td>
      <td>0.882318</td>
      <td>43</td>
      <td>0.170899</td>
      <td>6178</td>
      <td>6</td>
      <td>68</td>
      <td>22.5</td>
      <td>...</td>
      <td>-0.610256</td>
      <td>G1</td>
      <td>0.135907</td>
      <td>False</td>
      <td>1</td>
      <td>OPC</td>
      <td>6178</td>
      <td>305398</td>
      <td>13.918937</td>
      <td>6178</td>
    </tr>
    <tr>
      <th>CATAGACTCGATTATG-1_6178</th>
      <td>13373</td>
      <td>192511</td>
      <td>433</td>
      <td>0.224922</td>
      <td>66</td>
      <td>0.034284</td>
      <td>6178</td>
      <td>6</td>
      <td>68</td>
      <td>22.5</td>
      <td>...</td>
      <td>-3.405861</td>
      <td>G1</td>
      <td>0.063946</td>
      <td>False</td>
      <td>2</td>
      <td>DaN</td>
      <td>6178</td>
      <td>166685</td>
      <td>6.673376</td>
      <td>6178</td>
    </tr>
    <tr>
      <th>TCGTTTCCAATGCCTA-1_6178</th>
      <td>7726</td>
      <td>28196</td>
      <td>99</td>
      <td>0.351114</td>
      <td>44</td>
      <td>0.156051</td>
      <td>6178</td>
      <td>6</td>
      <td>68</td>
      <td>22.5</td>
      <td>...</td>
      <td>-0.613095</td>
      <td>G1</td>
      <td>0.078507</td>
      <td>False</td>
      <td>3</td>
      <td>Astro</td>
      <td>6178</td>
      <td>232259</td>
      <td>14.794874</td>
      <td>6178</td>
    </tr>
    <tr>
      <th>GTATGTTCAAGGACCA-1_6178</th>
      <td>7649</td>
      <td>28033</td>
      <td>108</td>
      <td>0.385260</td>
      <td>50</td>
      <td>0.178361</td>
      <td>6178</td>
      <td>6</td>
      <td>68</td>
      <td>22.5</td>
      <td>...</td>
      <td>-0.550458</td>
      <td>G1</td>
      <td>0.051028</td>
      <td>False</td>
      <td>3</td>
      <td>Astro</td>
      <td>6178</td>
      <td>163141</td>
      <td>16.353665</td>
      <td>6178</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>TTTGTGGCACGTTACA-1_5638</th>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5638</td>
      <td>9825</td>
      <td>11.113956</td>
      <td>5638</td>
    </tr>
    <tr>
      <th>TTTGTGGCATGAATAG-1_5638</th>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5638</td>
      <td>14257</td>
      <td>20.888272</td>
      <td>5638</td>
    </tr>
    <tr>
      <th>TTTGTGTTCAAGCGCC-1_5638</th>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5638</td>
      <td>20794</td>
      <td>5.329154</td>
      <td>5638</td>
    </tr>
    <tr>
      <th>TTTGTGTTCACGCGGT-1_5638</th>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5638</td>
      <td>1569</td>
      <td>21.080369</td>
      <td>5638</td>
    </tr>
    <tr>
      <th>TTTGTTGGTTGAGGTC-1_5638</th>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&lt;NA&gt;</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5638</td>
      <td>27298</td>
      <td>9.327705</td>
      <td>5638</td>
    </tr>
  </tbody>
</table>
<p>37281 rows × 24 columns</p>
</div>



We will save this object for later, after modeling latent spaces on the paired modality object.

### multiome_feature_selection

As with running scVI on the gene expression modality and spectral-MNN on chromatin accessibility, the most variable parameters should be used to infer UMAP and nearest neighbor latent spaces.

#### Gene expression layer


```python
# Select for the most variable genes
sc.pp.highly_variable_genes(
    rna, 
    n_top_genes=high_var_gene_num,
    subset=False,
    flavor = 'seurat_v3')
# Make a copy that only contains the highly variable genes
rna_hvg = rna[:, rna.var['highly_variable']].copy()
rna_hvg
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/legacy_api_wrap/__init__.py:88: UserWarning: `flavor='seurat_v3'` expects raw count data, but non-integers were found.
      return fn(*args_all, **kw)





    AnnData object with n_obs × n_vars = 23019 × 2000
        obs: 'n_genes_by_counts', 'total_counts', 'total_counts_rb', 'pct_counts_rb', 'total_counts_mt', 'pct_counts_mt', 'Sample_ID', 'Use_batch', 'Age', 'PMI', 'Sex', 'Brain_bank', 'Primary Diagnosis', 'S_score', 'G2M_score', 'phase', 'doublet_score', 'predicted_doublet', 'leiden', 'celltype'
        var: 'gene_ids', 'mt', 'rb', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
        uns: 'celltype_colors', 'hvg', 'leiden', 'leiden_colors', 'neighbors', 'umap'
        obsm: 'X_scvi', 'X_umap', 'padj_ulm', 'score_ulm'
        layers: 'counts', 'cpm', 'log-norm'
        obsp: 'connectivities', 'distances'




```python
# Verify the RNA layer is using counts (this is important for MultiVI)
if "counts" in rna_hvg.layers:
    rna_hvg.X = rna_hvg.layers["counts"].copy()
```

#### Chromatin accessibility layer


```python
# Select for the most variable peaks
snap.pp.select_features(atac, n_features=high_var_peak_num, inplace=True)
# Make a copy that only contains the highly variable peaks
atac_hvp = atac[:, atac.var["selected"].values].copy()
atac_hvp
```

    2026-11-03 18:42:19 - INFO - Selected 20000 features.





    AnnData object with n_obs × n_vars = 34934 × 20000
        obs: 'sample', 'n_fragment', 'tsse'
        var: 'count', 'selected'
        uns: 'spectral_eigenvalue'
        obsm: 'X_spectral', 'X_umap'
        obsp: 'distances'



Once we have the highly variable peaks and genes for both layers, we need to create 'dummy' nuclei that replace the nuclei that were filtered out in each modality. This will not impact the final multiome object latent space, but is required to run MultiVI.


```python
# Save modality specific barcodes
rna_barcodes = rna_hvg.obs_names
atac_barcodes = atac_hvp.obs_names

# Find the paired, RNA-, and ATAC-only barcodes
all_barcodes = rna_barcodes.union(atac_barcodes, sort=False)
paired_barcodes = rna_barcodes.intersection(atac_barcodes, sort=False)
rna_only_barcodes = rna_barcodes.difference(atac_barcodes, sort=False)
atac_only_barcodes = atac_barcodes.difference(rna_barcodes, sort=False)

print(f'{len(paired_barcodes)} paired nuclei')
print(f'{len(rna_only_barcodes)} nuclei with only gene expression')
print(f'{len(atac_only_barcodes)} nuclei with only chromatin accessibility')
```

    20672 paired nuclei
    2347 nuclei with only gene expression
    14262 nuclei with only chromatin accessibility


#### Fill in gaps in the individual modalities


```python
# Keep original RNA object; only create missing rows for ATAC-only nuclei
rna_missing = ad.AnnData(
    X=sp.csr_matrix((len(atac_only_barcodes), rna_hvg.n_vars), dtype=rna_hvg.X.dtype),
    obs=pd.DataFrame(index=atac_only_barcodes),
    var=rna_hvg.var.copy()
)
```


```python
# Keep original ATAC object; only create missing rows for RNA-only nuclei
atac_missing = ad.AnnData(
    X=sp.csr_matrix((len(rna_only_barcodes), atac_hvp.n_vars), dtype=atac_hvp.X.dtype),
    obs=pd.DataFrame(index=rna_only_barcodes),
    var=atac_hvp.var.copy()
)
```

#### MuData object ready for MultiVI


```python
# Create RNA and ATAC objects that have the same index, missing values
# will be filled the 'rna/atac_missing' objects
unified_rna = ad.concat([rna_hvg, rna_missing], merge="same")
unified_atac = ad.concat([atac_missing, atac_hvp], merge="same")

# Ensure only the nuclei with paired barcodes are included
unified_rna = unified_rna[all_barcodes, :]
unified_atac = unified_atac[all_barcodes, :]

# Create the object
hv_mdata = mu.MuData({"rna": unified_rna, "atac": unified_atac})
```

    /tmp/ipykernel_1641455/3460503119.py:3: UserWarning: Only some AnnData objects have `.raw` attribute, not concatenating `.raw` attributes.
      unified_rna = ad.concat([rna_hvg, rna_missing], merge="same")
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1471: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("var", axis=0, join_common=join_common)
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1322: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("obs", axis=1, join_common=join_common)



```python
# Metadata for plotting
paired_set = set(paired_barcodes)
rna_only_set = set(rna_only_barcodes)
atac_only_set = set(atac_only_barcodes)
hv_mdata.obs["modality_group"] = pd.Categorical(
    [
        "paired" if x in paired_set else
        "rna_only" if x in rna_only_set else
        "atac_only"
        for x in mdata.obs_names
    ],
    categories=["paired", "rna_only", "atac_only"]
)
```

#### Port over the `samples` column


```python
# Create sample key for batch-corrected MultiVI
hv_mdata.obs[sample_key] = [(x.split('-1_')[1]) for x in hv_mdata.obs_names]
# Update object
hv_mdata = hv_mdata.copy()
# Save object
hv_mdata
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1471: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("var", axis=0, join_common=join_common)
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1322: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("obs", axis=1, join_common=join_common)





<pre>MuData object with n_obs × n_vars = 37281 × 22000
  obs:	&#x27;modality_group&#x27;, &#x27;Sample_ID&#x27;
  2 modalities
    rna:	37281 × 2000
      var:	&#x27;gene_ids&#x27;, &#x27;mt&#x27;, &#x27;rb&#x27;, &#x27;highly_variable&#x27;, &#x27;highly_variable_rank&#x27;, &#x27;means&#x27;, &#x27;variances&#x27;, &#x27;variances_norm&#x27;
    atac:	37281 × 20000
      var:	&#x27;count&#x27;, &#x27;selected&#x27;</pre>



This result is an object with paired nuclei values for both modalities, where the are null values for nuclei in either modality that do not pass QC thresholds.

## Modeling

In earlier versions of scMAVERICS we used weighted-nearest neighbors to infer the multiome latent space. This unfortunately was incompatible with using unpaired modalities, which limited the number of nuclei in the final dataset. In the current version (version 1.2), [MultiVI](https://www.nature.com/articles/s41592-023-01909-9) version 1.4.2 is implemented. This allows for determining the latent space including all of the QC-passing data.

### multiome_model

Setting up the model parameters is important. There are methods to determine which are optimal for the input dataset, but in this case we will use parameters that have been shown to work for human CNS datasets.


```python
# Critical to specify the random number seed for reproducible stochasticity.
scvi.settings.seed = int(random_number_seed)
```

    INFO: [rank: 0] Seed set to 0
    2026-11-03 18:42:21 - INFO - [rank: 0] Seed set to 0



```python
# Setup SCVI on the data layer
scvi.model.MULTIVI.setup_mudata(
    hv_mdata,
    batch_key=sample_key,
    modalities={
        "rna_layer": "rna",
        "atac_layer": "atac"
    }
)

# Set up the model 
mvi_model = scvi.model.MULTIVI(
    hv_mdata,
    n_genes=hv_mdata.mod["rna"].n_vars,
    n_regions=hv_mdata.mod["atac"].n_vars,
)
```


```python
# Train the model
mvi_model.train(
    accelerator='cpu',
    max_epochs=max_epoch,
    lr=1e-3,
    early_stopping=True,
    batch_size=256
)
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/torch/nn/modules/linear.py:124: UserWarning: Initializing zero-element tensors is a no-op
      init.kaiming_uniform_(self.weight, a=math.sqrt(5))
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/lightning/fabric/plugins/environments/slurm.py:204: The `srun` command is available on your system but is not used. HINT: If your intention is to run Lightning on SLURM, prepend your python command with `srun` like so: srun python /data/catchingba/.conda/envs/multiome/lib/python3.11 ...
    INFO: GPU available: False, used: False
    2026-11-03 18:45:31 - INFO - GPU available: False, used: False
    INFO: TPU available: False, using: 0 TPU cores
    2026-11-03 18:45:31 - INFO - TPU available: False, using: 0 TPU cores
    INFO: 💡 Tip: For seamless cloud logging and experiment tracking, try installing [litlogger](https://pypi.org/project/litlogger/) to enable LitLogger, which logs metrics and artifacts automatically to the Lightning Experiments platform.
    2026-11-03 18:45:31 - INFO - 💡 Tip: For seamless cloud logging and experiment tracking, try installing [litlogger](https://pypi.org/project/litlogger/) to enable LitLogger, which logs metrics and artifacts automatically to the Lightning Experiments platform.
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/lightning/fabric/plugins/environments/slurm.py:204: The `srun` command is available on your system but is not used. HINT: If your intention is to run Lightning on SLURM, prepend your python command with `srun` like so: srun python /data/catchingba/.conda/envs/multiome/lib/python3.11 ...
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/lightning/pytorch/utilities/_pytree.py:21: `isinstance(treespec, LeafSpec)` is deprecated, use `isinstance(treespec, TreeSpec) and treespec.is_leaf()` instead.
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/lightning/pytorch/trainer/connectors/data_connector.py:434: The 'train_dataloader' does not have many workers which may be a bottleneck. Consider increasing the value of the `num_workers` argument` to `num_workers=63` in the `DataLoader` to improve performance.
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/lightning/pytorch/trainer/connectors/data_connector.py:434: The 'val_dataloader' does not have many workers which may be a bottleneck. Consider increasing the value of the `num_workers` argument` to `num_workers=63` in the `DataLoader` to improve performance.



    Training:   0%|          | 0/1000 [00:00<?, ?it/s]


    Monitored metric reconstruction_loss_validation did not improve in the last 50 records. Best score: 4787.161. Signaling Trainer to stop.



```python
elbo = mvi_model.history['elbo_train']
elbo['elbo_validation'] = mvi_model.history['elbo_validation']

plt.plot(elbo['elbo_train'])
```




    [<matplotlib.lines.Line2D at 0x7ffc70173610>]




    
![png](multiome_files/multiome_43_1.png)
    


### multiome_latent_transfer




```python
# Convert the cell barcode to the observable matrix X_scvi which neighbors and UMAP can be calculated from
mdata.obsm['X_multivi'] = mvi_model.get_latent_representation()
```


```python
mdata.obs["modality_group"] = hv_mdata.obs["modality_group"]
```


```python
# Calculate nearest neighbors and the UMAP from the X_scvi observable matrix
sc.pp.neighbors(mdata, use_rep='X_multivi')
sc.tl.umap(mdata, min_dist=0.3)
```


```python
# Calculate the leiden distance from the nearest neighbors, use a couple resolutions
sc.tl.leiden(mdata, key_added='leiden', flavor = 'igraph')
```


```python
# Transfer leiden clusters to rna and atac modalities
barcode2leiden = mdata.obs['leiden'].to_dict()
mdata.mod['rna'].obs['leiden'] = [barcode2leiden[x] for x in mdata.mod['rna'].obs_names]
mdata.mod['atac'].obs['leiden'] = [barcode2leiden[x] for x in mdata.mod['atac'].obs_names]
```


```python
mu.pl.umap(
    mdata,
    color = ['leiden'],
    frameon = False,
    legend_loc = 'on data'
)
```

    ... storing 'Sample_ID' as categorical



    
![png](multiome_files/multiome_51_1.png)
    



```python
mu.pl.umap(
    mdata,
    color = ['modality_group'],
    frameon = False,
)
```

    ... storing 'Sample_ID' as categorical



    
![png](multiome_files/multiome_52_1.png)
    


## Annotation

### multiome_cell_type_annotation

While the snRNAseq 


```python
# Create the DataFrame of canonical gene markers (This can be expanded)
marker_gene_df = pd.read_csv(work_dir + 'input/first_pass_genes.csv')
marker_gene_df
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
      <th>source</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Astro</td>
      <td>AQP4</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Astro</td>
      <td>GFAP</td>
    </tr>
    <tr>
      <th>2</th>
      <td>DaN</td>
      <td>SLC6A3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>DaN</td>
      <td>TH</td>
    </tr>
    <tr>
      <th>4</th>
      <td>EC</td>
      <td>CLDN5</td>
    </tr>
    <tr>
      <th>5</th>
      <td>EC</td>
      <td>PECAM1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>ExN</td>
      <td>SLC17A7</td>
    </tr>
    <tr>
      <th>7</th>
      <td>ExN</td>
      <td>SATB2</td>
    </tr>
    <tr>
      <th>8</th>
      <td>EpC</td>
      <td>TTR</td>
    </tr>
    <tr>
      <th>9</th>
      <td>EpC</td>
      <td>FOXJ1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>FB</td>
      <td>COL1A1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>FB</td>
      <td>LUM</td>
    </tr>
    <tr>
      <th>12</th>
      <td>InN</td>
      <td>GAD1</td>
    </tr>
    <tr>
      <th>13</th>
      <td>InN</td>
      <td>GAD2</td>
    </tr>
    <tr>
      <th>14</th>
      <td>MG</td>
      <td>ADAM28</td>
    </tr>
    <tr>
      <th>15</th>
      <td>MG</td>
      <td>APBB1IP</td>
    </tr>
    <tr>
      <th>16</th>
      <td>OPC</td>
      <td>PDGFRA</td>
    </tr>
    <tr>
      <th>17</th>
      <td>OPC</td>
      <td>VCAN</td>
    </tr>
    <tr>
      <th>18</th>
      <td>Oligo</td>
      <td>MOG</td>
    </tr>
    <tr>
      <th>19</th>
      <td>Oligo</td>
      <td>ST18</td>
    </tr>
    <tr>
      <th>20</th>
      <td>PC</td>
      <td>PDGFRB</td>
    </tr>
    <tr>
      <th>21</th>
      <td>PC</td>
      <td>KCNJ8</td>
    </tr>
    <tr>
      <th>22</th>
      <td>TC</td>
      <td>TRAC</td>
    </tr>
    <tr>
      <th>23</th>
      <td>TC</td>
      <td>CD3D</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Use only the rna modality
adata = mdata.mod['rna']

# Ensure the leiden clusters are the same between modalities
barcode2leiden = mdata.obs['leiden'].to_dict()
adata.obs['multiome_leiden'] = [barcode2leiden[x] for x in adata.obs_names]
adata
```




    AnnData object with n_obs × n_vars = 23019 × 36601
        obs: 'n_genes_by_counts', 'total_counts', 'total_counts_rb', 'pct_counts_rb', 'total_counts_mt', 'pct_counts_mt', 'Sample_ID', 'Use_batch', 'Age', 'PMI', 'Sex', 'Brain_bank', 'Primary Diagnosis', 'S_score', 'G2M_score', 'phase', 'doublet_score', 'predicted_doublet', 'leiden', 'celltype', 'multiome_leiden'
        var: 'gene_ids', 'mt', 'rb', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm'
        uns: 'celltype_colors', 'hvg', 'leiden', 'leiden_colors', 'neighbors', 'umap'
        obsm: 'X_scvi', 'X_umap', 'padj_ulm', 'score_ulm'
        layers: 'counts', 'cpm', 'log-norm'
        obsp: 'connectivities', 'distances'




```python
# Run over-represenation analysis based on cell markers
# provided in the marker_gene_df DataFrame.
dc.mt.ulm(
    adata,
    net = marker_gene_df,
    tmin=2,
    raw = False,
    layer = 'log-norm'
)
```


```python
# Convert the ORA AnnData object to numpy array to rank
score = dc.pp.get_obsm(adata, key="score_ulm")
df = dc.tl.rankby_group(adata=score, groupby='multiome_leiden', reference="rest", method="wilcoxon")
```


```python
# Apply the best ranked cell type to a cluster-celltype dictionary
annotation_dict = df.groupby('group').head(1).set_index('group')['name'].to_dict()

# Apply the dictionary to the MuData object
mdata.obs['celltype'] = [annotation_dict[clust] for clust in mdata.obs['leiden']]
```

    /tmp/ipykernel_1166891/890499784.py:2: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
      annotation_dict = df.groupby('group').head(1).set_index('group')['name'].to_dict()



```python
# Visualize the results
mu.pl.umap(
    mdata, 
    color = ['celltype', 'leiden'],
    legend_loc = 'on data',
    frameon = False
)
```

    ... storing 'celltype' as categorical



    
![png](multiome_files/multiome_61_1.png)
    



```python
# Port the cell types to the other layers
barcode2celltype = 
mdata.mod['atac'] = []
```


```python
sc.pl.dotplot(
    mdata.mod['rna'],
    marker_gene_df['target'],
    groupby = 'celltype',
    use_raw = False,
    layer = 'log-norm',
)
```

    ... storing 'multiome_leiden' as categorical



    
![png](multiome_files/multiome_63_1.png)
    


### multiome_pychromvar

One incredibly useful analysis we can run is the inference of transcription factor enrichment on a per nuclei basis. This is done using the package [pychromVAR](https://github.com/pinellolab/pychromVAR), a python port of [chromVAR](https://www.nature.com/articles/nmeth.4401) from the Greenleaf lab. We use the python port as the original packages were not built to scale to millions of nuclei/cells, and subsequently some of the methods are not very memory efficient. In the Snakemake pipeline this is addressed by chunking the dataset, but to finish in a reasonable amount of time (days) this still requires Tb of memory.

This is still an inferred modality, but it is much more reliable than using gene expression data as a proxy for transcription factor activity.


```python
# Get reference genome and sequences at peaks
pc.add_peak_seq(
    atac, 
    genome_file = ref_genome, 
    delimiter = ':|-' 
)
```

    100%|██████████| 608977/608977 [00:02<00:00, 265192.05it/s]



```python
# Correct for GC bias
pc.add_gc_bias(atac)
atac
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/pychromvar/preprocessing.py:117: FutureWarning: The method uns_keys is deprecated and will be removed in the future. Use uns instead of uns_keys. (e.g. `k in adata.uns` or `sorted(adata.uns)`)
      assert "peak_seq" in adata.uns_keys(), "Cannot find sequences, please first run add_peak_seq!"
    100%|██████████| 608977/608977 [00:02<00:00, 215346.86it/s]





    AnnData object with n_obs × n_vars = 34934 × 608977
        obs: 'sample', 'n_fragment', 'tsse'
        var: 'count', 'selected', 'gc_bias'
        uns: 'spectral_eigenvalue', 'peak_seq'
        obsm: 'X_spectral', 'X_umap'
        obsp: 'distances'




```python
reads_per_peak = np.log1p(atac.X.sum(axis=0)) / np.log(10)
reads_per_peak
```




    matrix([[1.75587484, 1.92427922, 2.521138  , ..., 2.50650497, 2.38201707,
             2.73798724]], shape=(1, 608977))




```python
atac.var['gc_bias'].values
```




    array([0.508, 0.55 , 0.466, ..., 0.482, 0.494, 0.6  ], shape=(608977,))




```python
"""pc.get_bg_peaks(mdata)"""
import scipy as sp
from pynndescent import NNDescent
niterations=50
n_jobs=-1
reads_per_peak = np.log1p(atac.X.sum(axis=0)) / np.log(10)

# here if reads_per_peak is a numpy matrix, convert it to array
if isinstance(reads_per_peak, np.matrix):
    reads_per_peak = np.squeeze(np.asarray(reads_per_peak))

mat = np.array([reads_per_peak, atac.var['gc_bias'].values])
chol_cov_mat = np.linalg.cholesky(np.cov(mat))
trans_norm_mat = sp.linalg.solve_triangular(
    a=chol_cov_mat, b=mat, lower=True).transpose()

index = NNDescent(trans_norm_mat, metric="euclidean",
                    n_neighbors=niterations, n_jobs=n_jobs)
knn_idx, _ = index.query(trans_norm_mat, niterations)

atac.varm['bg_peaks'] = knn_idx
```


```python
jdb_obj = jaspardb(release='JASPAR2026')
motifs = jdb_obj.fetch_motifs(collection='CORE', species=9606)
```


```python
pc.match_motif(atac, motifs=motifs)
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/pychromvar/match_motif.py:48: FutureWarning: The method uns_keys is deprecated and will be removed in the future. Use uns instead of uns_keys. (e.g. `k in adata.uns` or `sorted(adata.uns)`)
      assert "peak_seq" in adata.uns_keys(), "Cannot find sequences, please first run add_peak_seq!"
    100%|██████████| 608977/608977 [11:42<00:00, 866.82it/s]



```python
# Find how the individual nuclei TF motif enrichments deviate from the main score
chromvar = pc.compute_deviations(
    atac,
    n_jobs = 1
)
chromvar
```

    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/pychromvar/compute_deviations.py:40: FutureWarning: The method varm_keys is deprecated and will be removed in the future. Use varm instead of varm_keys. (e.g. `k in adata.varm` or `adata.varm.keys() | {'u'}`)
      assert "bg_peaks" in adata.varm_keys(
    2026-11-05 12:32:56 - INFO - computing expectation reads per cell and peak...
    2026-11-05 12:34:39 - INFO - computing observed motif deviations...
    2026-11-05 12:35:38 - INFO - computing background deviations...
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/legacy_api_wrap/__init__.py:88: FutureWarning: The dtype argument is deprecated and will be removed in late 2024.
      return fn(*args_all, **kw)





    AnnData object with n_obs × n_vars = 34934 × 857




```python
# Add the 
mdata.mod['chromvar'] = chromvar
```


```python
mdata.mod['chromvar'].var
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
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>MA0069.1.PAX6</th>
    </tr>
    <tr>
      <th>MA0071.1.RORA</th>
    </tr>
    <tr>
      <th>MA0074.1.RXRA::VDR</th>
    </tr>
    <tr>
      <th>MA0101.1.REL</th>
    </tr>
    <tr>
      <th>MA0107.1.RELA</th>
    </tr>
    <tr>
      <th>...</th>
    </tr>
    <tr>
      <th>MA2330.2.ZBTB24</th>
    </tr>
    <tr>
      <th>MA1589.3.ZNF140</th>
    </tr>
    <tr>
      <th>MA0088.3.ZNF143</th>
    </tr>
    <tr>
      <th>MA1976.3.ZNF320</th>
    </tr>
    <tr>
      <th>MA0116.2.ZNF423</th>
    </tr>
  </tbody>
</table>
<p>857 rows × 0 columns</p>
</div>




```python
mu.pl.umap(
    mdata,
    color = ['MA0473.4.ELF1', 'celltype'],
    cmap = 'PuOr_r',
    vcenter = 0
)
```

    ... storing 'celltype' as categorical



    
![png](multiome_files/multiome_76_1.png)
    


## Export the dataset


```python
mdata.write_h5mu(work_dir + 'data/08_multiome.h5mu', compression = 'gzip')
```

    ... storing 'celltype' as categorical
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1471: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("var", axis=0, join_common=join_common)
    /data/catchingba/.conda/envs/multiome/lib/python3.11/site-packages/mudata/_core/mudata.py:1322: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
      self._update_attr("obs", axis=1, join_common=join_common)

