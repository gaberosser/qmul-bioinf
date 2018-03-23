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
import numpy as np
import os
import pickle
from utils import output, setops
import multiprocessing as mp


def filter_cpm_by_group(cpm, groups, min_cpm=1):
    # filter
    over_min = (cpm > min_cpm).groupby(groups, axis=1).sum().astype(int)
    grp_size = the_groups.groupby(the_groups).size()
    return over_min.eq(grp_size).sum(axis=1) > 0



if __name__ == "__main__":
    # all samples to date (22-03-2018):
    pids = [
        '018', '019', '030', '031',
        '017', '050', '054', '061',
        '026', '052'
    ]
    subgroups = {
        'RTK I': ['018', '019', '030', '031'],
        'RTK II': ['017', '050', '054', '061'],
        'MES': ['026', '052']
    }

    # remove 061 due to limited methylation data (may be resolved shortly):
    # pids = [
    #     '018', '019', '030', '031',
    #     '017', '050', '054',
    #     '026', '052'
    # ]
    # subgroups = {
    #     'RTK I': ['018', '019', '030', '031'],
    #     'RTK II': ['017', '050', '054'],
    #     'MES': ['026', '052']
    # }

    # original 6 samples for comparison:
    # pids = [
    #     '018', '019', '031',
    #     '017', '050', '054',
    # ]
    # subgroups = {
    #     'RTK I': ['018', '019', '031'],
    #     'RTK II': ['017', '050', '054'],
    # }

    subgroup_set_colours = {
        'RTK I full': '#0d680f',
        'RTK II full': '#820505',
        'MES full': '#7900ad',
        'RTK I partial': '#6ecc70',
        'RTK II partial': '#d67373',
        'MES partial': '#cc88ea',
        'mixed': '#4C72B0',
        'specific': '#f4e842',
    }

    min_cpm = 1

    outdir = output.unique_output_dir("compare_de_gene_counts_s2", reuse_empty=True)

    # load our data

    obj = loader.load_by_patient(pids, include_control=True)

    # load reference data
    h9_obj = loader.load_references('GSE61794', tax_id=9606, source='star', strandedness='u')
    h1_obj = loader.load_references('GSE38993', tax_id=9606, source='star', strandedness='u', samples=['H1 NSC'])

    # combine
    obj = loader.MultipleBatchLoader([obj, h9_obj, h1_obj])

    # remove IPSC and rejected 061 samples for good
    idx = (
        (~obj.meta.index.str.contains('IPSC'))
        & (~obj.meta.index.isin(['DURA061_NSC_N1_P5', 'DURA061_NSC_N6_P4']))
    )
    obj.meta = obj.meta.loc[idx]
    obj.data = obj.data.loc[:, idx]
    obj.batch_id = obj.batch_id.loc[idx]

    refs = ['GIBCO', 'H1', 'H9']

    # we'll run everything with two different edgeR tests
    methods = ('GLM', 'QLGLM')

    res_1 = {}
    res_2 = {}
    step2_n_filter = {}
    res_3 = {}
    res_4 = {}
    step4_n_filter = {}

    pool = mp.Pool()
    jobs = {}

    # 1) Run GBM - iNSC DE separately for each, without filtering
    # this should give the same number of DE genes as we've seen before
    dat = obj.data.copy()

    for m in methods:

        this_res = {}
        for pid in pids:
            idx1 = (
                dat.columns.str.contains(pid)
            )
            the_data = dat.loc[:, idx1]
            the_groups = pd.Series('iNSC', index=the_data.columns)
            the_groups[the_groups.index.str.contains('GBM')] = 'GBM'
            the_comparison = ('GBM', 'iNSC')
            jobs[(m, pid, 'insc')] = pool.apply_async(
                differential_expression.run_one_de,
                args=(the_data, the_groups, the_comparison),
                kwds={'method': m}
            )

            for ref in refs:
                idx = (
                    (dat.columns.str.contains(pid) & dat.columns.str.contains('GBM'))
                    | dat.columns.str.contains(ref)
                )
                the_data = dat.loc[:, idx]
                the_groups = pd.Series(ref, index=the_data.columns)
                the_groups[the_groups.index.str.contains('GBM')] = 'GBM'
                the_comparison = ('GBM', ref)
                jobs[(m, pid, ref)] = pool.apply_async(
                    differential_expression.run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds={'method': m}
                )

    pool.close()
    pool.join()

    for (m, pid, typ), j in jobs.items():
        try:
            res_1.setdefault(m, {}).setdefault(pid, {})[typ] = j.get(1e4)
        except Exception as exc:
            print repr(exc)


    # 2) Run GBM - iNSC DE separately for each, filtering by CPM
    dat = obj.data.copy()

    pool = mp.Pool()
    jobs = {}

    for m in methods:

        for pid in pids:
            idx1 = (
                dat.columns.str.contains(pid)
            )
            the_data = dat.loc[:, idx1]
            the_groups = pd.Series('iNSC', index=the_data.columns)
            the_groups[the_groups.index.str.contains('GBM')] = 'GBM'

            keep = filter_cpm_by_group(
                the_data.divide(the_data.sum(), axis=1) * 1e6,
                the_groups,
                min_cpm=min_cpm
            )
            the_data = the_data.loc[keep]

            if m == methods[0]:
                step2_n_filter.setdefault(pid, {})['insc'] = (keep.size, keep.sum())

            the_comparison = ('GBM', 'iNSC')
            jobs[(m, pid, 'insc')] = pool.apply_async(
                differential_expression.run_one_de,
                args=(the_data, the_groups, the_comparison),
                kwds={'method': m}
            )

            for ref in refs:
                idx = (
                    (dat.columns.str.contains(pid) & dat.columns.str.contains('GBM'))
                    | dat.columns.str.contains(ref)
                )
                the_data = dat.loc[:, idx]
                the_groups = pd.Series(ref, index=the_data.columns)
                the_groups[the_groups.index.str.contains('GBM')] = 'GBM'

                keep = filter_cpm_by_group(
                    the_data.divide(the_data.sum(), axis=1) * 1e6,
                    the_groups,
                    min_cpm=min_cpm
                )
                the_data = the_data.loc[keep]

                if m == methods[0]:
                    step2_n_filter.setdefault(pid, {})[ref] = (keep.size, keep.sum())

                the_comparison = ('GBM', ref)
                jobs[(m, pid, ref)] = pool.apply_async(
                    differential_expression.run_one_de,
                    args=(the_data, the_groups, the_comparison),
                    kwds={'method': m}
                )

    pool.close()
    pool.join()

    for (m, pid, typ), j in jobs.items():
        try:
            res_2.setdefault(m, {}).setdefault(pid, {})[typ] = j.get(1e4)
        except Exception as exc:
            print repr(exc)


    # 3) Run GBM - iNSC DE using a pooled dispersion estimate, filtering first across all samples
    # AND
    # 4) As (3) but then adding additional filtering after the DE computation to remove genes that don't pass
    # the CPM requirement in the specific comparison being made
    dat = obj.data.copy()
    cpm = dat.divide(dat.sum(), axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) > 0
    dat = dat.loc[keep]
    print "Process 3. Pre-filtering all data together takes us from %d genes to %d (removing %d)." % (
        keep.size, keep.sum(), (~keep).sum()
    )

    # define groups now - they won't change (only the contrast will)

    the_groups = pd.Series(index=dat.columns)
    for pid in pids:
        the_groups[the_groups.index.str.contains('GBM') & the_groups.index.str.contains(pid)] = "GBM%s" % pid
        the_groups[the_groups.index.str.contains('NSC') & the_groups.index.str.contains(pid)] = "iNSC%s" % pid
    the_groups[the_groups.index.str.contains('GIBCO')] = "GIBCO"
    the_groups[the_groups.index.str.contains('H1')] = "H1"
    the_groups[the_groups.index.str.contains('H9')] = "H9"

    pool = mp.Pool()
    jobs = {}

    for m in methods:

        this_res1 = {}
        this_res2 = {}
        this_removed = {}

        for pid in pids:
            the_comparison = ('GBM%s' % pid, 'iNSC%s' % pid)
            jobs[(m, pid, 'insc')] = pool.apply_async(
                differential_expression.run_one_de,
                args=(dat, the_groups, the_comparison),
                kwds={'method': m}
            )

            for ref in refs:
                the_comparison = ('GBM%s' % pid, ref)
                jobs[(m, pid, ref.lower())] = pool.apply_async(
                    differential_expression.run_one_de,
                    args=(dat, the_groups, the_comparison),
                    kwds={'method': m}
                )

    pool.close()
    pool.join()

    for (m, pid, typ), j in jobs.items():
        try:
            res_3.setdefault(m, {}).setdefault(pid, {})[typ] = j.get(1e4)
        except Exception as exc:
            print repr(exc)

    # 4) filter again
    for m in methods:
        for pid in pids:
            idx1 = (
                dat.columns.str.contains(pid)
            )
            the_genes = res_3[m][pid]['insc'].index
            the_data = dat.loc[the_genes, idx1]
            the_groups = pd.Series('iNSC', index=the_data.columns)
            the_groups[the_groups.index.str.contains('GBM')] = 'GBM'

            keep = filter_cpm_by_group(
                the_data.divide(the_data.sum(), axis=1) * 1e6,
                the_groups,
                min_cpm=min_cpm
            )

            res_4.setdefault(m, {}).setdefault(pid, {})['insc'] = res_3[m][pid]['insc'].loc[keep]
            step4_n_filter.setdefault(m, {}).setdefault(pid, {})['insc'] = (keep.size, keep.sum())

            for ref in refs:
                idx = (
                    (dat.columns.str.contains(pid) & dat.columns.str.contains('GBM'))
                    | dat.columns.str.contains(ref)
                )
                the_genes = res_3[m][pid][ref].index
                the_data = dat.loc[the_genes, idx]
                the_groups = pd.Series(ref, index=the_data.columns)
                the_groups[the_groups.index.str.contains('GBM')] = 'GBM'

                keep = filter_cpm_by_group(
                    the_data.divide(the_data.sum(), axis=1) * 1e6,
                    the_groups,
                    min_cpm=min_cpm
                )

                res_4.setdefault(m, {}).setdefault(pid, {})[ref] = res_3[m][pid][ref].loc[keep]
                step4_n_filter.setdefault(m, {}).setdefault(pid, {})[ref] = (keep.size, keep.sum())




    # save the results - they don't take much space and make things much faster
    to_save = {
        'res_1': res_1,
        'res_2': res_2,
        'res_3': res_3,
        'res_4': res_4,
        'step2_n_filter': step2_n_filter,
        'step4_n_filter': step4_n_filter,
        'pids': pids,
        'subgroups': subgroups,
        'data': obj.data,
        'meta': obj.meta
    }
    fn_pkl = os.path.join(outdir, 'results.pkl')
    with open(fn_pkl, 'wb') as f:
        pickle.dump(to_save, f)
    print "Saved pickled results to %s" % fn_pkl

    # now let's look at the UpSet plot for each of these
    # first, we run with a reduced set of PIDs to match the previous analysis

    prev_pids = ['018', '019', '031', '017', '050', '054']
    prev_subgroups = {
        'RTK I': ['018', '019', '031'],
        'RTK II': ['017', '050', '054'],
    }
    sets_all = setops.full_partial_unique_other_sets_from_groups(prev_pids, prev_subgroups)

    set_colours = [
        ('RTK I full', {'sets': sets_all['full']['RTK I'], 'colour': '#0d680f'}),
        ('RTK I partial', {'sets': sets_all['partial']['RTK I'], 'colour': '#6ecc70'}),
        ('RTK II full', {'sets': sets_all['full']['RTK II'], 'colour': '#820505'}),
        ('RTK II partial', {'sets': sets_all['partial']['RTK II'], 'colour': '#d67373'}),
        ('Expanded core', {'sets': sets_all['mixed'], 'colour': '#4C72B0'}),
        ('Unique', {'sets': sets_all['specific'], 'colour': '#f4e842'})
    ]

    for m in methods:
        # UpsetR attribute plots
        data_for_upset1 = [res_1[m][pid].index for pid in prev_pids]  # this will be supplied to the function

        upset1 = venn.upset_set_size_plot(
            data_for_upset1,
            set_labels=prev_pids,
            set_colours=set_colours,
            min_size=10,
            n_plot=30,
            default_colour='gray'
        )
        upset1['figure'].savefig(os.path.join(outdir, "upset_%s_1_prev_pids.png" % m), dpi=200)

        data_for_upset2 = [res_2[m][pid].index for pid in prev_pids]  # this will be supplied to the function
        upset2 = venn.upset_set_size_plot(
            data_for_upset2,
            set_labels=prev_pids,
            set_colours=set_colours,
            min_size=10,
            n_plot=30,
            default_colour='gray'
        )
        upset2['figure'].savefig(os.path.join(outdir, "upset_%s_2_prev_pids.png" % m), dpi=200)

    # UpSet plots again, but this time including ALL the samples
    sets_all = setops.full_partial_unique_other_sets_from_groups(pids, subgroups)

    # create set colours
    # NB the order matters!

    set_colours = []
    for sg in subgroups:
        for x in ['full', 'partial']:
            k = "%s %s" % (sg, x)
            if sg in sets_all[x]:
                set_colours.append(
                    (k, {'sets': sets_all[x][sg], 'colour': subgroup_set_colours[k]})
                )
    set_colours.append(
        ('Expanded core', {'sets': sets_all['mixed'], 'colour': subgroup_set_colours['mixed']})
    )
    set_colours.append(
        ('Specific', {'sets': sets_all['specific'], 'colour': subgroup_set_colours['specific']})
    )

    for m in methods:
        data_for_upset1 = [res_1[m][pid].index for pid in pids]
        upset1 = venn.upset_set_size_plot(
            data_for_upset1,
            set_labels=pids,
            set_colours=set_colours,
            min_size=10,
            n_plot=30,
            default_colour='gray'
        )
        upset1['figure'].savefig(os.path.join(outdir, "upset_%s_1.png" % m), dpi=200)

        data_for_upset2 = [res_2[m][pid].index for pid in pids]  # this will be supplied to the function
        upset2 = venn.upset_set_size_plot(
            data_for_upset2,
            set_labels=pids,
            set_colours=set_colours,
            min_size=10,
            n_plot=30,
            default_colour='gray'
        )
        upset2['figure'].savefig(os.path.join(outdir, "upset_%s_2.png" % m), dpi=200)

        data_for_upset3 = [res_3[m][pid].index for pid in pids]
        upset3 = venn.upset_set_size_plot(
            data_for_upset3,
            set_labels=pids,
            set_colours=set_colours,
            min_size=10,
            n_plot=30,
            default_colour='gray'
        )
        upset3['figure'].savefig(os.path.join(outdir, "upset_%s_3.png" % m), dpi=200)

        data_for_upset4 = [res_4[m][pid].index for pid in pids]  # this will be supplied to the function
        upset4 = venn.upset_set_size_plot(
            data_for_upset4,
            set_labels=pids,
            set_colours=set_colours,
            min_size=10,
            n_plot=30,
            default_colour='gray'
        )
        upset4['figure'].savefig(os.path.join(outdir, "upset_%s_4.png" % m), dpi=200)