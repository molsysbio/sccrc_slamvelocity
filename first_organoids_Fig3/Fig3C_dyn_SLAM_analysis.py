'''
Main script to produce Fig 3 data.
'''

import numpy as np
import matplotlib.pyplot as pl
import scvelo as scv
import pandas as pd
import scanpy as sc

from scvelo_patient_utils import join_SLAM, find_roots

scv.settings.set_figure_params('scvelo')
scv.settings.verbosity = 4

data_path='/fast/scratch/users/peidlis_c/sodar_patient_organoid_data/'
signatures_path='/fast/work/users/peidlis_c/projects/sodar_patient_organoid_data/signatures/'

figures_path = '../figures/'
scv.settings.figdir=figures_path


if __name__ == "__main__":
    # Reads demuxed .h5 files from /demuxed, applies basic preprocessing, i.e.
    # filter by counts and mito read percentage, normalizes
    letters = ['C', 'E', 'W']
    donors = ['NCO', 'p009ot', 'p013ot']
    scv.settings.figdir = '/fast/work/users/peidlis_c/projects/sodar_patient_organoid_data/figures'

    panel_genes = ['FABP1', 'PHGR1', 'TFF3', 'MKI67', 'CD44']
    # my_best_genes=['EPCAM', 'KRT20', 'COL17A1', 'PHGR1', 'RPS3']


    ########## PREPARE DATA ############

    # UMAPS SLAM
    plot=True
    for donor in donors:
        adata = None

        for letter in letters:
            print('Slamming ', letter, donor)
            dat = scv.read(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'.h5')
            adata = join_SLAM(dat, letter)
            adata.write(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'SLAM.h5')
            if plot:
                scv.pl.scatter(adata, basis='umap', color=np.sum(adata.layers['unspliced'], axis=1), show=False, save='NB_AS_'+letter+'_'+donor+'_SLAMnewcounts.png', legend_loc='right margin')
                scv.pl.scatter(adata, basis='umap', color=np.sum(adata.layers['spliced'], axis=1), show=False, save='NB_AS_'+letter+'_'+donor+'_SLAMoldcounts.png', legend_loc='right margin')
                scv.pl.scatter(adata, basis='umap', color=panel_genes, show=False, save='NB_AS_'+letter+'_'+donor+'_SLAMmoments.png', legend_loc='right margin')

        # Aggregate conditions:
        print('Slamming Allconds ', donor)
        adata = None
        for letter in letters:
            if adata is None:
                adata = scv.read(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'SLAM.h5')
            else:
                bdata = scv.read(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'SLAM.h5')
                adata = adata.concatenate(bdata)
        # recompute stuff
        scv.pp.neighbors(adata)
        scv.pp.moments(adata)  # pool, normalize
        scv.tl.umap(adata)
        adata.write(data_path+'by_donors/NB_AS_'+'allconds_'+donor+'SLAM.h5')
        if plot:
            scv.pl.scatter(adata, basis='umap', color=np.sum(adata.layers['unspliced'], axis=1), show=False, save='NB_AS_allconds_'+donor+'_SLAMnewcounts.png', legend_loc='right margin')
            scv.pl.scatter(adata, basis='umap', color=np.sum(adata.layers['spliced'], axis=1), show=False, save='NB_AS_allconds_'+donor+'_SLAMoldcounts.png', legend_loc='right margin')
            scv.pl.scatter(adata, basis='umap', color=panel_genes, show=False, save='NB_AS_allconds_'+donor+'_SLAMmoments.png', legend_loc='right margin')

    ########## SLAM DYNAMICAL MODEL ############
    # dyn velo
    for donor in donors:
        # Split up by conditions and learn dynamics there
        for letter in letters:
            print(letter, donor)
            adata = scv.read(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'SLAMvelo.h5')

            roots = find_roots(adata, mode='plasticity', n=10)
            adata.uns['root_prior']=roots  # BUG, adata can not be saved with this

            dm = scv.tl.recover_dynamics(adata, var_names='velocity_genes', linear_action='check', fit_basal_transcription=False, root_prior=roots)
            scv.tl.velocity(adata, mode='dynamical', vkey='dyn_velo')
            scv.tl.velocity_graph(adata, vkey='dyn_velo')
            scv.tl.velocity_embedding(adata, basis='umap', vkey='dyn_velo')  # why not?
            scv.pl.velocity_embedding_stream(adata, color=panel_genes, vkey='dyn_velo', legend_loc='right margin', show=False, save='SLAMrooteddynvelopartial_NB_AS_'+letter+'_'+donor+'__panel.png')
            adata.write(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'SLAMvelo.h5')



        # Allconds
        # adata = scv.read(data_path+'by_donors/NB_AS_'+'allconds_'+donor+'SLAMvelo.h5')
        # roots = find_roots(adata, mode='plasticity', n=10)
        # adata.uns['root_prior']=roots  # BUG, adata can not be saved with this
        # dm = scv.tl.recover_dynamics(adata, var_names='velocity_genes', linear_action='check', fit_basal_transcription=False, root_prior=roots)
        # scv.tl.velocity(adata, mode='dynamical', vkey='dyn_velo')
        # scv.tl.velocity_graph(adata, vkey='dyn_velo')
        # scv.tl.velocity_embedding(adata, basis='umap', vkey='dyn_velo')
        # scv.pl.velocity_embedding_stream(adata, color='condition', vkey='dyn_velo', legend_loc='right margin', show=False, save='SLAMrooteddynvelo_NB_AS_allconds_'+donor+'_condition.png')
        # scv.pl.velocity_embedding_stream(adata, color=panel_genes, vkey='dyn_velo', legend_loc='right margin', show=False, save='SLAMrooteddynvelo_NB_AS_allconds_'+donor+'_panel.png')
        # for letter in letters:
        #     print(letter, donor)
        #     bdata = adata
        #     scv.pl.velocity_embedding_stream(adata, vkey='dyn_velo', legend_loc='right margin', show=False, save='SLAMrooteddynvelo_NB_AS_allconds_'+donor+'.png')
        # #adata.write(data_path+'by_donors/NB_AS_'+'allconds_'+donor+'SLAMvelo.h5')

    print('done')
