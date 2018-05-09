from rnaseq import loader
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


if __name__ == "__main__":
    pids = ['017', '018', '019', '030', '031', '044', '050', '054']
    bin_edges = np.linspace(0, 20000, 101)

    obj = loader.load_by_patient(pids, include_control=False)
    dat = obj.data.loc[:, obj.data.columns.str.contains('GBM')]
    dat = dat.loc[dat.sum(axis=1) > 0]
    mean_dat = dat.mean(axis=1)

    grouped_by_avg = []
    m = []
    for i in range(len(bin_edges) - 1):
        idx = (mean_dat > bin_edges[i]) & (mean_dat <= bin_edges[i + 1])
        m.append(dat.loc[idx].mean(axis=1).mean())
        if idx.sum() > 2:
            grouped_by_avg.append(dat.loc[idx].var(axis=1))
        else:
            grouped_by_avg.append(None)
        # grouped_by_avg.append(dat.loc[(dat > bin_edges[i]) & (dat <= bin_edges[i + 1])].values)

