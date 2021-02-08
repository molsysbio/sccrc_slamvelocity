'''
Reads the HTOdemux and cellSNP+vireo results and applies them to the looms.

Input:
- HTOdemux result: data_path+'NB_AS_'+letter+'/demuxed/'+letter+'_donor_ids_seurat.csv'
- vireo result: data_path+'NB_AS_'+letter+'/demuxed/'+letter+'_donor_ids_SNP.tsv'
- velocyto looms: data_path+'NB_AS_'+letter+'/20200328_velocyto/NB_AS_'+letter+'.loom'

Output:
- demuxed data: data_path+'NB_AS_'+letter+'/demuxed/NB_AS_'+letter+'_demuxed.h5'
'''

import numpy as np
import matplotlib.pyplot as pl
import scvelo as scv
import pandas as pd
import scanpy as sc

scv.settings.set_figure_params('scvelo')
scv.settings.verbosity = 4

data_path='/fast/scratch/users/peidlis_c/sodar_patient_organoid_data/'

def apply_demux(adata, letter, data_path='./data/'):
    annot = {'C': {
              'Hashtag10': 'p009ot',
              'Hashtag11': 'p013ot',
              'Hashtag12': 'NCO',
              },
        'E': {
              'Hashtag7': 'p009ot',
              'Hashtag8': 'p013ot',
              'Hashtag9': 'NCO',
              },
        'W': {
              'Hashtag4': 'p009ot',
              'Hashtag5': 'p013ot',
              'Hashtag6': 'NCO',
              }
        }

    SNPdemux = pd.read_csv(data_path+'NB_AS_'+letter+'/demuxed/'+letter+'_donor_ids_SNP.tsv', sep='\t')
    SNPdemux = SNPdemux.set_index('cell', drop=True)
    new_index = [index[:-2] for index in SNPdemux.index]  # clean index names
    SNPdemux.index=new_index

    adata.obs['SNPdemux'] = SNPdemux['donor_id']

    seuratdemux = pd.read_csv(data_path+'NB_AS_'+letter+'/demuxed/'+letter+'_donor_ids_seurat.csv', sep=',')
    for key in annot[letter].keys():
        seuratdemux[seuratdemux==key]=annot[letter][key]

    for name in pd.unique(seuratdemux['x']):
        if name==name:
            if '_' in name:
                seuratdemux[seuratdemux==name]='doublet'

    adata.obs['seuratdemux'] = seuratdemux
    adata.obs['seuratdemux'][pd.isna(adata.obs['seuratdemux'])]='Negative'
    return adata

def identify(adata, letter):
    annot = {'C': {
              'Hashtag10': 'p009ot',
              'Hashtag11': 'p013ot',
              'Hashtag12': 'NCO',
              },
        'E': {
              'Hashtag7': 'p009ot',
              'Hashtag8': 'p013ot',
              'Hashtag9': 'NCO',
              },
        'W': {
              'Hashtag4': 'p009ot',
              'Hashtag5': 'p013ot',
              'Hashtag6': 'NCO',
              }
        }
    # apply seurat id to SNP
    donor_names = np.array(list(annot[letter].values()))
    donor_ids = np.array(['donor0', 'donor1', 'donor2'])
    mat=np.zeros((len(donor_ids), len(donor_names)))
    for i, sn in enumerate(donor_ids):
        for j, se in enumerate(donor_names):
            mat[i][j] = len(set(np.where(adata.obs['SNPdemux']==sn)[0]).intersection(set(np.where(adata.obs['seuratdemux']==se)[0])))
    # pl.imshow(mat)
    donors_dict = dict(list(zip(donor_ids[np.argmax(mat, axis=0)], donor_names)))
    for key in donors_dict.keys():
        adata.obs['SNPdemux'][adata.obs['SNPdemux']==key]=donors_dict[key]
    return adata

if __name__ == "__main__":
    # Takes the original loom files from SODAR, demuxes them using the seurat
    # /BP label (_donor_ids_seurat.csv) and SNP information (_donor_ids_SNP.tsv)
    # saves result as _demuxed.h5 containing demux info in obs
    letters = ['C', 'E', 'W']
    adatas = [scv.read(data_path+'NB_AS_'+letter+'/20200328_velocyto/NB_AS_'+letter+'.loom') for letter in letters]
    for cdata, letter in zip(adatas, letters):
        new_index = [index[25:-1] for index in cdata.obs.index]  # clean index names
        cdata.obs_names=new_index
        apply_demux(cdata, letter, data_path=data_path)
        identify(cdata, letter)
        cdata.write(data_path+'NB_AS_'+letter+'/demuxed/NB_AS_'+letter+'_demuxed.h5')
