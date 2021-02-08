'''
This script splits the demuxed libraries into donors, preprocessed them and
subsequently splits each into a single .h5.

Input: data_path+'NB_AS_'+letter+'/demuxed/NB_AS_'+letter+'_demuxed.h5'
Output: data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'.h5'
'''

import numpy as np
import matplotlib.pyplot as pl
import scvelo as scv
import pandas as pd
import scanpy as sc

scv.settings.set_figure_params('scvelo')
scv.settings.verbosity = 4

data_path='/fast/scratch/users/peidlis_c/sodar_patient_organoid_data/'
signatures_path='/fast/work/users/peidlis_c/sodar_patient_organoid_data/signatures/'

colors = {
    'p009ot': "#E69F00",
    'p013ot': "#0072B2",
    'G1':  "#94b6d2",
    'S': "#dc4040",
    'G2M' : "#7aa77f",
    'Stem' :"#FEC98DFF",
    'Cycling TA' :"#F1605DFF",
    'Secretory TA' :"#CD4071FF",
    'Immature Enterocytes 1':"#7AD151FF",
    'Immature Enterocytes 2':"#6BCC3EFF",
    'Enterocyte Progenitors':'#00f5b4',
    'Enterocytes' :"#009E73",
    'Enteroendocrine': '#451077FF',
    'Immature Goblet' :"#9F2F7FFF",
    'Goblet':"#451077FF",
    'M cells' :"#56B4E9",
    'Tuft':"#FDE725FF",
    'NA' : "grey",
    'Cycling B' :"#FEC98DFF",
    'B cells' :"#FEC98DFF",
    'Plasma' :"#F1605DFF",
    'Follicular' :"#9F2F7FFF",
    'GCs' :"#451077FF",
    'Treg' :"#0072B2",
    'Tconv' :"#56B4E9",
    'CD4+ T':"#D55E00",
    'CD8+ T' :"#CC79A7",
    'Mast' :"#D55E00",
    'DC1' :"#009E73",
    'DC2' :"#F0E442",
    'ILCs' :"#0072B2",
    'Monocytes' :"#F0E442",
    'Macroph.' :"#009E73",
    'Crypt FBs' :"#FEC98DFF",
    'Villus FBs' :"#F1605DFF",
    'Infl. FBs' :"#56B4E9",
    'Myofib.' :"#009E73",
    'Endothelial':"#F0E442",
    'PC venules' :"#0072B2",
    'Pericytes' :"#D55E00",
    'Glia' :"#CC79A7",
    'normal' : "steelblue",
    'tumor' : "red",
    'PDO' : "steelblue"
}



def pp(adata, min_counts=1000, max_perc_mito=0.25, louvain_res=1.8, reg_cc=False):
    # Mito reads
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['raw_counts'] = adata.X.sum(axis=1).A1
    adata = adata[adata.obs.percent_mito < max_perc_mito, :].copy()

    # simple pp
    sc.pp.filter_cells(adata, min_counts=min_counts)  # Nils sagt Daten (organoide) sind spitze, also 1000 min counts
    sc.pp.filter_genes(adata, min_cells=2)
    sc.pp.normalize_total(adata)

    # cc score
    cell_cycle_genes = [x.strip() for x in open(signatures_path+'cell_cycle_genes/regev_lab_cell_cycle_genes.txt')]
    # Split into 2 lists
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    if reg_cc:
        sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
        sc.pp.scale(adata)

    # annotate plasticity (~ #genes expressed)
    activity = np.array(adata.X > np.mean(adata.X, axis=0))
    plasticity=np.sum(activity, axis=1)
    adata.obs['plasticity'] = plasticity

    # precomps
    scv.pp.neighbors(adata)
    scv.pp.moments(adata)
    scv.tl.umap(adata)

    # clustering
    scv.tl.louvain(adata, resolution=louvain_res)
    return adata

if __name__ == "__main__":
    # Reads demuxed .h5 files from /demuxed, applies basic preprocessing, i.e.
    # filter by counts and mito read percentage, normalizes

    # Florians curated colon epithelial cell type signature, derived from regev colon markers
    tab = pd.read_csv(signatures_path + 'cell_type_markers/_data__tab_diff_markers_epi.tsv', sep='\t')
    tab = tab[tab.p_val_adj<0.01]
    tab = tab[np.abs(tab.avg_logFC)> 0.75]
    signatures = dict([(ct, tab.gene.values[tab.cell_type_epi_custom==ct]) for ct in pd.unique(tab.cell_type_epi_custom)])

    letters = ['C', 'E', 'W']
    donors = ['NCO', 'p009ot', 'p013ot']
    scv.settings.figdir = '/fast/work/users/peidlis_c/sodar_patient_organoid_data/figures'
    # Note: Never analyse NCO and patient tumors together.

    panel_genes = ['FABP1', 'PHGR1', 'TFF3', 'MKI67', 'CD44']

    # prep and save data
    recalc = True
    if recalc:
        # Split up conditions:
        for letter in letters:
            cdata = scv.read(data_path+'NB_AS_'+letter+'/demuxed/NB_AS_'+letter+'_demuxed.h5')
            for donor in donors:
                print(letter, donor)
                adata = cdata[cdata.obs['SNPdemux']==donor].copy()
                adata.var_names_make_unique()
                adata=pp(adata, min_counts=5000, max_perc_mito=.25)
                # ct_annotate(adata, signatures, rescore=True)
                adata.write(data_path+'NB_AS_'+letter+'/processed/NB_AS_'+letter+'_'+donor+'.h5')

        # Aggregate conditions:
        for donor in donors:
            adata = None
            for letter in letters:
                dat = scv.read(data_path+'NB_AS_'+letter+'/demuxed/NB_AS_'+letter+'_demuxed.h5')
                dat = dat[dat.obs['SNPdemux']==donor].copy()  # select donor
                dat.obs['condition']=letter
                # aggregate
                if adata is None:
                    adata=dat
                else:
                    adata=adata.concatenate(dat, index_unique=None)
            adata.var_names_make_unique()
            adata=pp(adata, min_counts=5000, max_perc_mito=.25)
            # ct_annotate(adata, signatures, rescore=True)
            adata.write(data_path+'by_donors/NB_AS_'+'allconds_'+donor+'.h5')

    print('done')
