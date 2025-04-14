# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import scanpy as sc
import pycistarget
from pycisTopic.lda_models import evaluate_models, un_cgs_models_mallet

with open(snakemake.input.merged_cistopic_object, "rb") as f:
    cistopic_obj = pickle.load(f)

os.environ['MALLET_MEMORY'] = '200G'
# Configure path Mallet
mallet_path="/data/CARD_singlecell/SN_atlas/data/pycisTopic/Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=48,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="/data/CARD_singlecell/SN_atlas/data/pycisTopic/mallet",
    save_path="/data/CARD_singlecell/SN_atlas/data/pycisTopic/mallet",
    mallet_path=mallet_path,
)

cistopic_obj.add_LDA_model(model)

pickle.dump(
    models,
    open(snakemake.output.cistopic_model, "wb")
)

pickle.dump(
    cistopic_obj,
    open(snakemake.output.merged_cistopic_object, "wb")
)