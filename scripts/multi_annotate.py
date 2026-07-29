import muon as mu

mdata = mu.read(work_dir + '/atlas/06_polished_multiome_rna.h5mu')
del mdata.mod['rna'].obs['atlas_identifier']

# 
leiden_dict = mdata.obs['leiden'].to_dict()
mdata.mod['rna'].obs['mdata_leiden'] = [leiden_dict[x] for x in mdata.mod['rna'].obs_names]