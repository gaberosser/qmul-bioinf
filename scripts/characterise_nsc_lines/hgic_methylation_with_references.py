from methylation import loader, process
from plotting import clustering
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations
import os
from settings import LOCAL_DATA_DIR



def hc_plot_dendrogram_vary_n_gene(
        data,
        row_colours,
        mad = None,
        n_ftr=(1000, 2000, 3000, 5000, 10000),
        metric='correlation'
):
    """
    For each value in n_gene_arr, plot a dendrogram showing the result of hierarchical clustering of the data using
    that many genes (selected in descending MAD order)
    :param data: Cols are samples, rows are genes (or similar)
    :param row_colours: As passed to dendrogram routine
    :param n_gene_arr: The values to test
    :return:
    """
    if mad is None:
        mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
    fig_dict = {}
    for ng in n_ftr:
        the_dat = data.loc[mad.index[:ng]]
        d = clustering.dendrogram_with_colours(
            the_dat,
            row_colours,
            fig_kws={'figsize': (5.5, 10)},
            vertical=False,
            metric=metric
        )
        fig_dict[ng] = d
    return fig_dict



if __name__ == "__main__":
    norm_method = 'bmiq'
    outdir = unique_output_dir("methylation_insc_characterisation", reuse_empty=True)
    # load 12 patients iNSC, 4 iPSC
    pids = ['017', '018', '019', '030', '031', '026', '044', '049', '050', '052', '054', '061']

    patient_obj = loader.load_by_patient(pids, norm_method=norm_method)

    refs = [
        ('Kim et al.', loader.gse38216(norm_method=norm_method, samples=['H9 ESC 1', 'H9 ESC 2', 'H9 NPC 1', 'H9 NPC 2'])),
        ('Morey et al.', loader.gse67283(norm_method=norm_method, samples=['NPC_wt'])),
        ('Zimmerlin et al.', loader.gse65214(norm_method=norm_method)),
        ('Encode EPIC', loader.encode_epic(norm_method=norm_method)),
        # ('Encode 450k', loader.encode_450k(norm_method=norm_method)),
    ]

    to_aggr = [
        (r'H9 ESC [12]', 'H9 ESC'),
        (r'H9 NPC [12]', 'H9 NPC'),
    ]

    # TODO make this nicer
    r = refs[0][1]
    for srch, repl in to_aggr:
        idx = r.data.columns.str.contains(srch)
        new_col = r.data.loc[:, idx].mean(axis=1)
        r.data = r.data.loc[:, ~idx]
        r.data.insert(r.data.shape[1], repl, new_col)
        new_meta_row = r.meta.loc[idx].iloc[0]
        new_meta_row.name = repl
        r.meta = r.meta.loc[~idx]
        r.meta = r.meta.append(new_meta_row)

    # rename based on study
    for nm, t in refs:
        new_idx = ["%s (%s)" % (i, nm) for i in t.meta.index]
        t.meta.index = new_idx
        t.data.columns = new_idx

    ref_obj = loader.loader.MultipleBatchLoader([t[1] for t in refs])

    obj = loader.loader.MultipleBatchLoader([patient_obj, ref_obj])

    bdat = obj.data
    mdat = process.m_from_beta(bdat)
    # bmad = transformations.median_absolute_deviation(bdat).sort_values(ascending=False)
    mmad = transformations.median_absolute_deviation(mdat).sort_values(ascending=False)

    row_colours_all = pd.DataFrame('gray', index=mdat.columns, columns=[''])
    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = '#96daff'
    row_colours_all.loc[row_colours_all.index.str.contains(r'GBM')] = '#fff89e'
    row_colours_all.loc[row_colours_all.index.str.contains(r'NS27Z')] = '#fff89e'

    row_colours_all.loc[row_colours_all.index.str.contains('H9')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('H7')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('H1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('ESO3')] = '#ff7777'

    row_colours_all.loc[row_colours_all.index.str.contains('neuron')] = '#ccebc5'
    row_colours_all.loc[row_colours_all.index.str.contains(r'[Aa]strocyte')] = '#e78ac3'

    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = '#fdc086'  # orange

    n_ftr = [1000, 2000, 3000, 5000, 10000, 20000, 50000, 100000, 1000000]
    plt_dict = hc_plot_dendrogram_vary_n_gene(mdat, row_colours_all, mad=mmad, n_ftr=n_ftr)
    for ng, x in plt_dict.items():
        fname = "dendrogram_M_corr_top%d_by_mad.{ext}" % ng
        x['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=200)
        # x['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)