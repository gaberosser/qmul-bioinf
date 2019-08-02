from load_data import rnaseq_data
import os
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np


if __name__ == '__main__':
    star_count_dir = rnaseq_data.wtchg_p170582.params['star']['count_dirs'][0]
    cuff_count_dir = rnaseq_data.wtchg_p170582.params['star']['cufflinks']['count_dirs'][0]
    meta_fn = rnaseq_data.wtchg_p170582.meta_files[0]

    obj_star = rnaseq_data.StarCountLoader(count_dir=star_count_dir, meta_fn=meta_fn, strandedness='r')
    obj_cuff = rnaseq_data.CufflinksGeneLoader(count_dir=cuff_count_dir, meta_fn=meta_fn)

    star_fpkm = obj_star.get_fpkm()
    cuff_fpkm = obj_cuff.data

    common_index= star_fpkm.index.intersection(cuff_fpkm.index)
    star_fpkm = star_fpkm.loc[common_index]
    cuff_fpkm = cuff_fpkm.loc[common_index]

    samples = [u'DURA031_NSC_N44_P3', u'DURA030_NSC_N9_P2', u'DURA019_NSC_N5C1_P2'][::-1]
    for sn in samples:
        x_star = star_fpkm.loc[:, sn]
        x_cuff = cuff_fpkm.loc[:, sn]
        lr = stats.linregress(np.log(x_star + 1), np.log(x_cuff + 1))
        print "Sample %s, R^2 = %.3f" % (sn, lr.rvalue ** 2)
        fig = plt.figure(num=sn)
        plt.scatter(
            np.log(x_star + 1),
            np.log(x_cuff + 1),
        )
        plt.title(sn)
        plt.xlabel("My FPKM estimate (log)")
        plt.ylabel("Cufflinks FPKM (log)")
