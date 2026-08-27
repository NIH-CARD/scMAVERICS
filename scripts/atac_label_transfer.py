rule atac_label_transfer:
    input:
        merged_rna_anndata = work_dir+'/atlas/07_polished_anndata_rna.h5ad',
        merged_atac_anndata = work_dir+'/atlas/03_filtered_anndata_atac.h5ad'
    output:
        merged_atac_anndata = work_dir+'/atlas/04_annot_anndata_atac.h5ad'
    params:
        pseudobulk_param = 'celltype'
    script:
        'scripts/atac_label_transfer.py'
