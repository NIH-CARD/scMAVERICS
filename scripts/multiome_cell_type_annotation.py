import numpy as np
import pandas as pd

import scanpy as sc
import decoupler as dc
import muon as mu

mdata = mu.read_h5mu(snakemake.input.multiome_object)

