"""
21-03-2018
Purpose of this script: compute DE gene counts under a number of different conditions.
1) Pooled dispersion estimation vs separate
2) Sample choice (in case of pooled dispersion estimation)
3) Statistical test
"""


from rnaseq import loader, general, differential_expression
from matplotlib import pyplot as plt
import seaborn as sns
from plotting import venn
import pandas as pd


if __name__ == "__main__":
    pids = [
        '018', '019', '030', '031',
        '017', '050', '054', '061',
        '026', '052'
    ]
    min_cpm = 1

    obj = loader.load_by_patient(pids, include_control=True)

    # remove IPSC and rejected 061 samples for good
    idx = (~obj.meta.index.str.contains('IPSC')) & (~obj.meta.index.isin(['DURA061_NSC_N1_P5', 'DURA061_NSC_N6_P4']))
    obj.meta = obj.meta.loc[idx]
    obj.data = obj.data.loc[:, idx]
    obj.batch_id = obj.batch_id.loc[idx]

    # we'll run everything with two different edgeR tests

    methods = ('GLM', 'QLGLM')

    res_1 = {}
    res_2 = {}
    step2_n_filter = {}
    res_3 = {}
    res_4 = {}
    step4_n_filter = {}

    for m in methods:

        # 1) Run GBM - iNSC DE separately for each, without filtering
        # this should give the same number of DE genes as we've seen before
        dat = obj.data.copy()
        idx = (dat.columns.str.contains('GBM') | dat.columns.str.contains('NSC')) & (~dat.columns.str.contains('GIBCO'))
        dat = dat.loc[:, idx]

        this_res = {}
        for pid in pids:
            the_data = dat.loc[:, dat.columns.str.contains(pid)]
            the_groups = pd.Series('iNSC', index=the_data.columns)
            the_groups[the_groups.index.str.contains('GBM')] = 'GBM'
            the_comparison = ('GBM', 'iNSC')
            this_res[pid] = differential_expression.run_one_de(
                the_data,
                the_groups,
                the_comparison,
                method=m
            )
        res_1[m] = this_res

        # 2) Run GBM - iNSC DE separately for each, filtering by CPM
        dat = obj.data.copy()
        idx = (dat.columns.str.contains('GBM') | dat.columns.str.contains('NSC')) & (~dat.columns.str.contains('GIBCO'))
        dat = dat.loc[:, idx]

        this_res = {}
        for pid in pids:
            the_data = dat.loc[:, dat.columns.str.contains(pid)]
            the_groups = pd.Series('iNSC', index=the_data.columns)
            the_groups[the_groups.index.str.contains('GBM')] = 'GBM'

            # filter
            cpm = the_data.divide(the_data.sum(), axis=1) * 1e6
            over_min = (cpm > min_cpm).groupby(the_groups, axis=1).sum().astype(int)
            grp_size = the_groups.groupby(the_groups).size()
            keep = over_min.eq(grp_size).sum(axis=1) > 0
            the_data = the_data.loc[keep]
            if m == methods[0]:
                step2_n_filter[pid] = (keep.size, keep.sum())

            the_comparison = ('GBM', 'iNSC')
            this_res[pid] = differential_expression.run_one_de(
                the_data,
                the_groups,
                the_comparison,
                method=m
            )
        res_2[m] = this_res

        # 3) Run GBM - iNSC DE using a pooled dispersion estimate, filtering first across all samples
        # AND
        # 4) As (3) but then adding additional filtering after the DE computation to remove genes that don't pass
        # the CPM requirement in the specific comparison being made
        dat = obj.data.copy()
        idx = (dat.columns.str.contains('GBM') | dat.columns.str.contains('NSC')) & (~dat.columns.str.contains('GIBCO'))
        dat = dat.loc[:, idx]
        cpm = dat.divide(dat.sum(), axis=1) * 1e6
        keep = (cpm > min_cpm).sum(axis=1) > 0
        dat = dat.loc[keep]
        if m == methods[0]:
            print "Process 3. Pre-filtering all data together takes us from %d genes to %d (removing %d)." % (
                keep.size, keep.sum(), (~keep).sum()
            )

        the_groups = pd.Series(index=dat.columns)
        for pid in pids:
            the_groups[the_groups.index.str.contains('GBM') & the_groups.index.str.contains(pid)] = "GBM%s" % pid
            the_groups[the_groups.index.str.contains('NSC') & the_groups.index.str.contains(pid)] = "iNSC%s" % pid

        this_res1 = {}
        this_res2 = {}
        this_removed = {}
        for pid in pids:
            the_comparison = ('GBM%s' % pid, 'iNSC%s' % pid)
            this_res1[pid] = differential_expression.run_one_de(
                dat,
                the_groups,
                the_comparison,
                method=m
            )

            # 4) filter again
            the_data = dat.loc[this_res1[pid].index, dat.columns.str.contains(pid)]
            grp = pd.Series('iNSC', index=the_data.columns)
            grp[grp.index.str.contains('GBM')] = 'GBM'

            cpm = the_data.divide(the_data.sum(), axis=1) * 1e6
            over_min = (cpm > min_cpm).groupby(grp, axis=1).sum().astype(int)
            grp_size = grp.groupby(grp).size()
            keep = over_min.eq(grp_size).sum(axis=1) > 0

            this_res2[pid] = this_res1[pid].loc[keep]
            this_removed[pid] = (keep.size, keep.sum())

        res_3[m] = this_res1
        res_4[m] = this_res2
        step4_n_filter[m] = this_removed
