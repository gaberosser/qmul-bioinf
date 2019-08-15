from plotting import bar, common, pie
from methylation import loader, dmr, process
import pandas as pd
from statsmodels.sandbox.stats import multicomp
from utils import output, setops, genomics, log
import multiprocessing as mp
import os
import collections
import pickle
import numpy as np
from scipy import stats, cluster
import matplotlib
from matplotlib import pyplot as plt, patches
from matplotlib.colors import Normalize
from matplotlib import cm
from sklearn.neighbors import KernelDensity

import seaborn as sns
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts

from settings import HGIC_LOCAL_DIR, LOCAL_DATA_DIR, INTERMEDIATE_DIR
logger = log.get_console_logger()

# FIXME: this is a hack to try to avoid non-thread safe TKinter issue, which results in a segfault when we try
# to use multiprocessing. This requires ipython NOT to be run with --matplotlib
matplotlib.use('Agg')


if __name__ == "__main__":
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    dmr_params['n_jobs'] = mp.cpu_count()

    # set this to True if output bed files are required (this is quite slow due to the large number of combinations)
    write_bed_files = False

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')

    me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES_GIC)
    me_data = me_obj.data
    me_meta = me_obj.meta
    me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    delta_m = {}
    for col1 in me_meta.index:
        delta_m[col1] = {}
        for col2 in me_meta.index:
            if col1 == col2:
                continue
            delta_m[col1][col2] = me_data[col2] - me_data[col1]

    # median values
    dm_median = pd.DataFrame(index=me_meta.index.copy(), columns=me_meta.index.copy(), dtype=float)
    dm_median.index.name = 'A'
    dm_median.columns.name = 'B'
    for col1 in me_meta.index:
        for col2 in me_meta.index:
            if col1 == col2:
                continue
            dm_median.loc[col1, col2] = delta_m[col1][col2].median()


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.tick_top()
    sns.heatmap(dm_median, ax=ax)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "median_delta_m_pairwise.png"), dpi=200)

    # # array of KDEs
    # edges = np.linspace(-5, 5, 40)
    # fig, axs = plt.subplots(nrows=len(me_meta.index), ncols=len(me_meta.index), sharex=True, sharey=True)
    # for i, col1 in enumerate(me_meta.index):
    #     for j, col2 in enumerate(me_meta.index):
    #         if col1 == col2:
    #             axs[i, j].set_visible(False)
    #             continue
    #         axs[i, j].hist(delta_m[col1][col2], edges)
    #         axs[i, j].set_xlim([-5, 5])
    #         if i == 0:
    #             axs[i, j].set_xlabel(col2, rotation=90)
    #             axs[i, j].xaxis.set_label_position('top')
    #         if j == 0:
    #             axs[i, j].set_ylabel(col1, rotation=0)
    #             axs[i, j].yaxis.set_ticks([])
    #
