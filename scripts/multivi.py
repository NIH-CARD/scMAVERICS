import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import snapatac2 as snap
import scvi
import scipy.sparse as sp

# Load multiome dataset
hvg_mdata = mu.read(sys.argv[1])
# Setup
scvi.model.MULTIVI.setup_mudata(
    hvg_mdata,
    batch_key = sys.argv[2],
    modalities={
        "rna_layer": "rna",
        "atac_layer": "atac"
    }
)

# Train
mvi_model = scvi.model.MULTIVI(
    hvg_mdata,
    n_genes=hvg_mdata.mod["rna"].n_vars,
    n_regions=hvg_mdata.mod["atac"].n_vars,
)

mvi_model.train(
    accelerator="cpu",
    max_epochs=1000,
    lr=1e-3,
    early_stopping=True,
    batch_size=256
)

# Latent space transfer
hvg_mdata.obsm["X_multivi"] = mvi_model.get_latent_representation()
# Compute nearest neighbors
sc.pp.neighbors(hvg_mdata, use_rep="X_multivi")

# Compute leiden
sc.tl.leiden(hvg_mdata, key_added='leiden')

# Compute umap
sc.tl.umap(hvg_mdata)

hvg_mdata.write(sys.argv[4], compression='gzip')
mvi_model.save(sys.argv[5], overwrite=True)
