import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy as sc
import matplotlib.backends.backend_pdf
import scvelo as scv
from tqdm import tqdm, trange

scv.settings.set_figure_params('scvelo', dpi_save=100, dpi=80, transparent=True)
scv.settings.verbosity = 0
sc.settings.verbosity = 0

data_path='/fast/scratch/users/peidlis_c/sodar_patient_organoid_data/'
signatures_path='/fast/work/users/peidlis_c/projects/sodar_patient_organoid_data/signatures/'

letters = ['C', 'E', 'W']
donors = ['NCO', 'p009ot', 'p013ot']
scv.settings.figdir = '/fast/work/users/peidlis_c/sodar_patient_organoid_data/figures'

scoring_method='seurat'
new_umaps=False

# What to plot
Regev=False
CRC_stems=True
Hallmarks=False
custom_sigs=False
speed=False
metadata=False

def plot_signature(signature, name):
    from matplotlib import cm as cm
    my_cmap = cm.get_cmap('magma')
    for i in range(20):
        my_cmap.colors[i]=[.7]*3

    pdf = matplotlib.backends.backend_pdf.PdfPages(
          './figures/NB_AS_'+name+'.pdf')
    for ids in tqdm(signature.keys(), leave=None):
        fig, ax = pl.subplots(3, 3, figsize=[14, 12])

        genes = signature[ids]
        genes = [genes] if np.isscalar(genes) else genes

        # get values
        for j, donor in enumerate(donors):
            for i, letter in enumerate(letters):
                adata = adatas[i, j]
                if scoring_method == 'seurat':
                    sc.tl.score_genes(adata, genes, score_name=ids, use_raw=False)
                else:
                    hits = np.where(np.isin(adata.var_names, genes))[0]
                    if len(hits) > 0:
                        c = adata.layers['Ms'][:, hits]
                        cp=np.mean(c,axis=1)
                        adata.obs[ids]=cp

        # normalize
        vmin=0
        vmax=1
        for j, donor in enumerate(donors):
            for i, letter in enumerate(letters):
                adata = adatas[i, j]
                vmax = max(vmax, np.max(adata.obs[ids]))
                vmin = min(vmin, np.min(adata.obs[ids]))

        # plot
        for j, donor in enumerate(donors):
            for i, letter in enumerate(letters):
                adata = adatas[i, j]
                scv.pl.scatter(adata, color=ids, show=False, color_map=my_cmap,
                                   ax=ax[i, j], title=ids + ' ' + donor + ' ' + letter, vmin=vmin, vmax=vmax, colorbar= j==2 & i==2)

        pdf.savefig(fig)
    pdf.close()

if __name__ == "__main__":

    # load data
    print('loading data...')
    adatas = np.empty((3,3), dtype=object)
    for j, donor in enumerate(tqdm(donors)):
        for i, letter in enumerate(letters):
            adata = scv.read(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'SLAM.h5')
            if new_umaps:
                bdata = adata.copy()
                sc.pp.filter_genes_dispersion(bdata, n_top_genes=2000)
                sc.tl.pca(bdata, n_comps=50)
                sc.pp.neighbors(bdata)
                sc.tl.umap(bdata)
                adata.obsm = bdata.obsm
            adatas[i,j] = adata
    print('finished loading.')

    # # Regev
    if Regev:
    	print('Processing regev...')
        name='regev_signature'
        tab = pd.read_excel(signatures_path+'cell_type_markers/regev-colon-mmc2.xlsx', index_col=1, header=0)
        tab = tab[tab.pvalD<0.01]
        tab = tab[np.abs(tab.pvalC)> 0.75]
        signature=dict([(ct, tab.index.values[tab['ident']==ct]) for ct in pd.unique(tab['ident'])])
        plot_signature(signature, name)

    # CRC stem signatures
    if CRC_stems:
	print('Processing CRC stem signatures...')
	tab = pd.read_excel(signatures_path+'cell_type_markers/CRC-related_stem_cell_signatures.xlsx', index_col=None, header=0, skiprows=[1])
	name = 'CRC-related_stem_cell_signatures_check_old'
	signature = dict([(ids, tab[ids][~pd.isna(tab[ids])].values) for ids in tab.columns])
	plot_signature(signature, name)

    # # Hallmark signatures
    if Hallmarks:
	print('Processing hallmarks...')
	tab = pd.read_csv(signatures_path+'hallmarks/h.all.v6.2.symbols.gmt', sep='\t', index_col=None, header=None).T.drop(1)
	tab.columns=tab.iloc[0]
	tab = tab.drop(0)
	signature = dict([(ids[9:], tab[ids][~pd.isna(tab[ids])].values) for ids in tab.columns])
	name = 'hallmark_annotation'
	plot_signature(signature, name)

    # Florian custom signatures
    if custom_sigs:
	print('Processing epi_all signature...')
	tab = pd.read_csv(signatures_path + 'cell_type_markers/cell_type_markers_epi.tsv', sep='\t')
	tab = tab[tab.p_val_adj<0.01]
	tab = tab[np.abs(tab.avg_logFC)> 0.75]
	signature=dict([(ct, tab.gene.values[tab.cell_type_epi==ct]) for ct in pd.unique(tab.cell_type_epi)])
	plot_signature(signature, name='colonoid_cancer_uhlitz_markers_epi_thelistwithOLFM4')

    # SPEED pathways
    if speed:
	print('Processing progeny signature...')
	tab = pd.read_csv(signatures_path+'targets/progeny_pathways.csv', sep='\t', header=0)
	signature = dict([(ids, tab[tab['#Pathway']==ids].Gene_Name.values) for ids in pd.unique(tab['#Pathway'])])
	plot_signature(signature, name='progeny_pathways')

    # metadata plots
    if metadata:
	print('Processing metadata...')
	name = 'metadata'
	pdf = matplotlib.backends.backend_pdf.PdfPages(
	  './figures/NB_AS_'+name+'.pdf')
	for j, donor in enumerate(donors):
	for i, letter in enumerate(letters):
	    adata = adatas[i,j]
	    scv.pl.scatter(adata, basis='umap', color=['percent_mito', 'raw_counts', 'phase', 'plasticity'],
	                   show=False)#, save='NB_AS_'+letter+'_'+donor+'_UMAP.png', legend_loc='right margin')
	    pdf.savefig(pl.gcf())

	adata = adatas_all[j]
	scv.pl.scatter(adata, basis='umap', color=['percent_mito', 'condition', 'raw_counts', 'phase', 'plasticity'],
	               show=False)#, save='NB_AS_allconds_'+donor+'_UMAP.png', legend_loc='right margin')
	pdf.savefig(pl.gcf())
	pdf.close()
