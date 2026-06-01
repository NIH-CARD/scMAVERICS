import numpy as np
import pandas as pd
import scanpy as sc
import dreampy as dp

# Import pseudobulked data as an AnnData object
# This needs to be Sample-Celltype by gene
pb = sc.read_h5ad(snakemake.input.pseudo_rna)

# Parameters to slice 
celltype_param = snakemake.params.celltype_params
cell_grouping = snakemake.params.celltypes
diagnosis_param = snakemake.params.diagnosis_param
diagnosis_control = snakemake.params.diagnosis_values

# Change parameter names from the Decoupler hard-coded values to
# Dreampy hard-coded values
dc_pb.obs['n_cells'] = dc_pb.obs['psbulk_cells']
dc_pb.obs['assays'] = dc_pb.obs[celltype_param]

# Preprocessing
pb = dp.filter_samples(pb, min_cells=10, min_samples=3)
pb = dp.compute_tmm_factors(pb, assay_col=celltype_param)
assays = dp.filter_by_expr(pb, assay_col=celltype_param)

# These variables are needed in the Dreampy analysis
nf = pb.obs["norm_factors"].values
geo_mean = np.exp(np.mean(np.log(nf)))

# Assign the split cell types to each 
assays_dict = dp.filter_by_expr(pb, assay_col=celltype_param)

# Create a list of all possible comparisons
comparison_combinations = []
for i, condition_1 in enumerate(diagnosis_control):
    for j, condition_2 in enumerate(diagnosis_control):
        if i > j:
            comparison_combinations.append([condition_1, condition_2])

for celltype in assays_dict.keys():
    print(f'\n-- {celltype} --')
    assay_pb = assays_dict[celltype]
    for comparison in comparison_combinations:
        
        if 'control' in comparison:
            disease_name = comparison[0]
        else:
            disease_name = f'{comparison[0]} vs. {comparison[1]}'
        print(f'\n-- {disease_name} --')
        comparison_assay_pb = assay_pb[assay_pb.obs[diagnosis_param].isin(comparison)]

        # log2CPM
        comparison_assay_pb = dp.log2cpm(comparison_assay_pb)

        if comparison_assay_pb.shape[0] != 0:

            # Voom weights
            comparison_assay_pb = dp.estimate_weights(
                comparison_assay_pb, formula=snakemake.params.formula, n_jobs=-1,
            )

            # Mixed-model fitting
            fit = dp.fit_models(
                comparison_assay_pb,
                formula=snakemake.params.formula,
                reml=True,
                n_jobs=snakemake.threads,
                assay_name=celltype,
            )
            print(f"  fit: {fit}")

            # eBayes
            fit_eb = dp.ebayes(fit)
            print(f"  eBayes: s2_prior={fit_eb.s2_prior:.4f}, "
                f"df_prior={fit_eb.df_prior:.1f}")

            # Results table
            results = dp.get_results(fit_eb, assay_name=celltype)
            results[celltype_param] = celltype
            results[diagnosis_param] = disease_name

            all_results = pd.concat([all_results, results])
all_results.to_csv(snakemake.output.output_DGE_data, compression='gzip')