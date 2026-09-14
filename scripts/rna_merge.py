import anndata as ad
import pandas as pd
import scanpy as sc
import scipy

# Create dictionary of sample name to file, then filter with relevant diagnosis
sample_loc = dict(zip(snakemake.params.samples, snakemake.input.rna_anndata))
adatas = {key: sc.read_h5ad(sample_loc[key]) for key in snakemake.params.samples}

# Concatenate with names of sample
adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas=adatas
    )

# Produce a sparce counts layer
adata.layers['counts'] = scipy.sparse.csr_matrix(adata.layers['counts'].copy())
adata.layers['cpm'] = scipy.sparse.csr_matrix(adata.layers['cpm'].copy())
adata.layers['log-norm'] = scipy.sparse.csr_matrix(adata.layers['log-norm'].copy())

## Estimating sex and outputting table of XIST vs. RPS4Y1 expression
# Pseudobulk data
pdata = dc.pp.pseudobulk(
    adata,
    sample_col=snakemake.params.samples,
    layer='counts',
    mode='sum'
)

# Dropping unreliable pseudobulk samples
dc.pp.filter_samples(pdata, min_cells=snakemake.params.min_cells)

# Normalizing 
pdata.layers["counts"] = pdata.X.copy()
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)

# Getting expression of XIST and RPS4Y1
xist_idx = pdata.var_names.get_loc('XIST')
rps4y1_idx = pdata.var_names.get_loc('RPS4Y1')

xist_expr = pdata.X[:, xist_idx]
rps4y1_expr = pdata.X[:, rps4y1_idx]
xist_expr = xist_expr.toarray().flatten() if hasattr(xist_expr, 'toarray') else np.asarray(xist_expr).flatten()
rps4y1_expr = rps4y1_expr.toarray().flatten() if hasattr(rps4y1_expr, 'toarray') else np.asarray(rps4y1_expr).flatten()

q1_xist, q3_xist = np.percentile(xist_expr, [25, 75])
q1_rps4y1, q3_rps4y1 = np.percentile(rps4y1_expr, [25, 75])

xist_high = xist_expr > q3_xist
xist_low = xist_expr < q1_xist
rps4y1_high = rps4y1_expr > q3_rps4y1
rps4y1_low = rps4y1_expr < q1_rps4y1

# Predict sex
sex_pred = np.array(['Unknown'] * pdata.n_obs, dtype=object)
sex_pred[xist_high & rps4y1_low] = 'Female'
sex_pred[xist_low & rps4y1_high] = 'Male'

pdata.obs['predicted_sex'] = pd.Categorical(sex_pred, categories=['Female', 'Male', 'Unknown'])
pdata.obs['XIST_expr'] = xist_expr
pdata.obs['RPS4Y1_expr'] = rps4y1_expr

sex_map = pdata.obs.set_index(snakemake.params.samples)['predicted_sex']
adata.obs['predicted_sex'] = adata.obs[snakemake.params.samples].map(sex_map)

# pull reported sex from adata.obs onto the pseudobulk samples for coloring
reported_sex_map = adata.obs.drop_duplicates(subset=snakemake.params.samples).set_index(snakemake.params.samples)[snakemake.params.sex]
pdata.obs[snakemake.params.sex] = pdata.obs[snakemake.params.samples].map(reported_sex_map)

# Make a plot
color_map = {'Male': 'tab:blue', 'Female': 'tab:red'}
colors = pdata.obs[snakemake.params.sex].map(color_map).fillna('gray')

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(pdata.obs['XIST_expr'], pdata.obs['RPS4Y1_expr'], c=colors, s=60, edgecolor='black', linewidth=0.5)

ax.axvline(q1_xist, color='gray', linestyle='--', linewidth=0.8)
ax.axvline(q3_xist, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(q1_rps4y1, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(q3_rps4y1, color='gray', linestyle='--', linewidth=0.8)

unknown_mask = pdata.obs['predicted_sex'] == 'Unknown'
for sample_id, x, y in zip(
    pdata.obs.loc[unknown_mask, snakemake.params.samples],
    pdata.obs.loc[unknown_mask, 'XIST_expr'],
    pdata.obs.loc[unknown_mask, 'RPS4Y1_expr']
):
    ax.annotate(sample_id, (x, y), fontsize=7, xytext=(4, 4), textcoords='offset points')

handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:red', markersize=8, label='Female (reported)'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:blue', markersize=8, label='Male (reported)'),
]
if colors.eq('gray').any():
    handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='No Sex listed in metadata'))
ax.legend(handles=handles, loc='best', fontsize=8)

ax.set_xlabel('XIST normalized expression (log1p CPM)')
ax.set_ylabel('RPS4Y1 normalized expression (log1p CPM)')
ax.set_title('Pseudobulk sex marker expression by sample')

fig.tight_layout()
fig.savefig(snakemake.output.sex_prediction_plot, format='pdf')
plt.close(fig)

# Write out the unfiltered dataset
adata.write_h5ad(
    filename=snakemake.output.merged_rna_anndata, 
    compression='gzip'
    )