import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.weightstats import ztest as ztest

# Open the AnnData object
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Input parameters
sample_key = snakemake.params.sample_key
disease_param = snakemake.params.disease_param
seperating_cluster = snakemake.params.separating_cluster

# Check all the samples that they have more than one cell type in them
sample_total = adata.obs.groupby(['Sample']).count()['total_counts'].to_dict()

#age_sample = adata.obs[[sample_key, 'Age']].drop_duplicates()#.reindex(adata.obs['Sample_ID'].drop_duplicates().to_list())
#age_sample.index = adata.obs[sample_key].drop_duplicates().to_list()
#age_sample = age_sample['Age'].to_dict()

sample_disease = adata.obs[['Sample', disease_param]]
sample_disease.index = sample_disease['Sample']
sample2disease = sample_disease[disease_param].to_dict()

sample_num_cells = adata.obs.groupby([sample_key, seperating_cluster]).count().reset_index()[[sample_key, seperating_cluster, 'total_counts']]
sample_num_cells.columns = [sample_key, seperating_cluster, 'cell type counts']
#sample_num_cells['Age'] = [age_sample[x] for x in sample_num_cells[sample_key].to_list()]
sample_num_cells['total cells'] = [sample_total[x] for x in sample_num_cells[sample_key].to_list()]
sample_num_cells['cell type counts'] = sample_num_cells['cell type counts']
sample_num_cells[disease_param] = [sample2disease[x] for x in sample_num_cells[sample_key].to_list()]
sample_num_cells['fraction of sample'] = sample_num_cells['cell type counts'] / sample_num_cells['total cells']
sample_num_cells['log10 fraction of sample'] = np.log10(sample_num_cells['fraction of sample'])

plt.figure(figsize=(6, 2))
g = sns.boxplot(
        data=sample_num_cells.sort_values(by=disease_param),
        x=seperating_cluster,
        y='log10 fraction of sample',
        hue=disease_param,
        hue_order=list(snakemake.params.separating_value_dict.keys()),
        palette = list(snakemake.params.separating_value_dict.values(),)
        boxprops={'edgecolor': 'none'},
        flierprops={"markersize": 1, 'marker': '.'},
        showcaps=False
)
plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
plt.title('Sample cell type composition by primary diagnosis', fontweight='bold', fontsize=5)
plt.ylabel('percent of sample')
plt.xlabel('')
plt.tick_params(length=0)
g.spines[['right', 'top']].set_visible(False)
g.spines[['bottom', 'left']].set_edgecolor('#cccccc')
plt.yticks([0, -1, -2, -3, -4], [100, 10, 1, 0.1, 0.01])
plt.tight_layout(pad=.4)
plt.savefig(snakemake.output.fraction_boxplot)

# Store the stats in an array
cell_stats = []

for cell in cell_types:
    for disease in snakemake.params.diseases:
        sample_num_cell_type = sample_num_cells[(sample_num_cells[seperating_cluster] == cell)]
        control_values = sample_num_cell_type[sample_num_cell_type[disease_param] == snakemake.params.control]['fraction of sample']
        disease_values = sample_num_cell_type[sample_num_cell_type[disease_param] == disease]['fraction of sample']
        cell_stats.append([cell, disease, ztest(control_values, disease_values)[0], ztest(control_values, disease_values)[1]])
        #print(cell, disease, ztest(disease_values, control_values))
cell_stat_df = pd.DataFrame(cell_stats)
cell_stat_df.columns = [seperating_cluster, disease_param, 'z-test stat', 'z-test p-value']
cell_stat_df.to_csv(snakemake.output.corrected_ztest_results, index=False)
cell_stat_df
