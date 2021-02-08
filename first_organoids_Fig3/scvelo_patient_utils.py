'''
Util functions for SLAM_analysis.py.
'''

import numpy as np
import matplotlib.pyplot as pl
import scvelo as scv
import pandas as pd
import scanpy as sc


scv.settings.set_figure_params('scvelo')
scv.settings.verbosity = 4

data_path='/fast/scratch/users/peidlis_c/sodar_patient_organoid_data/'
signatures_path='/fast/work/users/peidlis_c/projects/sodar_patient_organoid_data/signatures/'


def find_roots(adata, mode='plasticity', n=10):
    if mode == 'stem_cluster':
        # Use n random stem cluster cells as prior if available
        if 'Stem' in pd.unique(adata.obs['annot_test']):
            roots = np.random.choice(np.where(adata.obs['annot_test']=='Stem')[0], n)
        else:
            raise Exception('no root found')
            roots=None
    elif mode == 'signature':
        # Use cells that are top n at Stem signature
        ctms = pd.read_csv(signatures_path+'cell_type_markers/cell_type_markers_epi.tsv', sep='\t', index_col=1, header=0)
        key='cell_type_epi'
        cell_type = 'Stem'
        marker=ctms.index[ctms[key]==cell_type]
        hits=np.where(np.isin(adata.var_names, marker))[0]
        c = adata.layers['Ms'][:, hits]
        cp=np.mean(c,axis=1)
        roots = np.argsort(cp)[-n:]
    elif mode == 'plasticity':
        # select the n cells with highest plasticity (~ #genes expressed)
        if 'plasticity' not in adata.obs.keys():
            activity = np.array(adata.X > np.mean(adata.X, axis=0))
            plasticity=np.sum(activity, axis=1)
            adata.obs['plasticity'] = plasticity
        roots = np.argsort(adata.obs['plasticity'])[-n:]
    # adata.uns['root_prior']=roots
    return roots


def join_SLAM(bdata, letter):
    for time in ['old', 'new']:
        adata= sc.read_10x_mtx('/fast/scratch/groups/ag_bluethgen/count'+letter+'/'+time+'_matrix', make_unique=False)
        adata.var_names_make_unique()

        # intersect vars
        bdata = bdata[:, np.isin(bdata.var_names, adata.var_names)].copy()
        adata = adata[:, np.isin(adata.var_names, bdata.var_names)].copy()

        # intersect obs (filtered)
        new_index = [index[:-2] for index in adata.obs_names]  # clean index names
        adata.obs_names=new_index
        adata = adata[np.isin(adata.obs_names, bdata.obs_names), :].copy()

        # align
        adata = adata[:, np.argsort(adata.var_names)].copy()
        bdata = bdata[:, np.argsort(bdata.var_names)].copy()

        adata = adata[np.argsort(adata.obs_names), :].copy()
        bdata = bdata[np.argsort(bdata.obs_names), :].copy()

        print('adata ', adata.n_obs, adata.n_obs)
        print('bdata ', bdata.n_obs, bdata.n_obs)
        # add
        bdata.layers[time]=adata.X

    # overwrite u and s, then redo pool, then normalize, but keep umap. No filtering
    bdata.layers['unspliced'] = bdata.layers['new']
    bdata.layers['spliced'] = bdata.layers['old']
    scv.pp.neighbors(bdata)
    scv.pp.moments(bdata)  # pool, normalize
    return bdata
