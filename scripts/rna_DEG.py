import numpy as np
import pandas as pd
import scanpy as sc
import dreampy as dp
import matplotlib.pyplot as plt
import seaborn as sns

pb = sc.read_h5ad('/data/CARD_singlecell/SN_atlas/atlas/pseudobulked_rna.h5ad')

celltype_param = snakemake.params.celltype_params
cell_grouping = snakemake.params.celltypes
diagnosis_param = snakemake.params.diagnosis_param
diagnosis_control = snakemake.params.diagnosis_values

pb = dp.aggregate_pseudobulk(adata, groupby=[cell_grouping, "Sample_ID"])

# Preprocessing
pb = dp.filter_samples(pb, min_cells=10, min_samples=3)
pb = dp.compute_tmm_factors(pb, assay_col=cell_grouping)
assays = dp.filter_by_expr(pb, assay_col=cell_grouping)

nf = pb.obs["norm_factors"].values
geo_mean = np.exp(np.mean(np.log(nf)))

assays_dict = dp.filter_by_expr(pb, assay_col=cell_grouping)



comparison_combinations = []
for i, condition_1 in enumerate(diagnosis_control):
    for j, condition_2 in enumerate(diagnosis_control):
        if i > j:
            comparison_combinations.append([condition_1, condition_2])

FORMULA = "~ Primary Diagnosis + Age + Sex + (1|psbulk_counts) + (1|psbulk_cells) + (1|Use_batch) + (1|Brain_bank)"

for comparison in comparison_combinations:
        print(f'\n-- {comparison[0]} vs. {comparison[1]} --')
        if comparison[0] == 'control':
            disease_name = comparison[1]
        else:
            disease_name = f'{comparison[0]} vs. {comparison[1]}'
        print(f'\n-- {disease_name} --')
        comparison_assay_pb = assay_pb[assay_pb.obs[diagnosis_param].isin(comparison)]

        # log2CPM
        comparison_assay_pb = dp.log2cpm(comparison_assay_pb)

        if comparison_assay_pb.shape[0] != 0:

            # Voom weights
            comparison_assay_pb = dp.estimate_weights(
                comparison_assay_pb, formula=FORMULA, n_jobs=-1,
            )

            #voom_plot(assay_pb, assay_name)
            #plt.show()

            # Mixed-model fitting
            fit = dp.fit_models(
                comparison_assay_pb,
                formula=FORMULA,
                reml=True,
                n_jobs=-1,
                assay_name=assay_name,
            )
            print(f"  fit: {fit}")

            # eBayes
            fit_eb = dp.ebayes(fit)
            print(f"  eBayes: s2_prior={fit_eb.s2_prior:.4f}, "
                f"df_prior={fit_eb.df_prior:.1f}")

            # Results table
            results = dp.get_results(fit_eb, assay_name=assay_name)
            results[celltype_param] = assay_name
            results[diagnosis_param] = disease_name
            all_results = pd.concat([all_results, results])