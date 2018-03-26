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
    grp_size = groups.groupby(groups).size()
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
    # h1_obj = loader.load_references('GSE38993', tax_id=9606, source='star', strandedness='u', samples=['H1 NSC'])

    # combine
    obj = loader.MultipleBatchLoader([obj, h9_obj])

    # remove IPSC and rejected 061 samples for good
    idx = (
        (~obj.meta.index.str.contains('IPSC'))
        & (~obj.meta.index.isin(['DURA061_NSC_N1_P5', 'DURA061_NSC_N6_P4']))
    )
    obj.meta = obj.meta.loc[idx]
    obj.data = obj.data.loc[:, idx]
    obj.batch_id = obj.batch_id.loc[idx]

    refs = ['GIBCO', 'H9']

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
            jobs[(m, pid, 'iNSC')] = pool.apply_async(
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
                step2_n_filter.setdefault(pid, {})['iNSC'] = (keep.size, keep.sum())

            the_comparison = ('GBM', 'iNSC')
            jobs[(m, pid, 'iNSC')] = pool.apply_async(
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
            jobs[(m, pid, 'iNSC')] = pool.apply_async(
                differential_expression.run_one_de,
                args=(dat, the_groups, the_comparison),
                kwds={'method': m}
            )

            for ref in refs:
                the_comparison = ('GBM%s' % pid, ref)
                jobs[(m, pid, ref)] = pool.apply_async(
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
            the_genes = res_3[m][pid]['iNSC'].index
            the_data = dat.loc[the_genes, idx1]
            the_groups = pd.Series('iNSC', index=the_data.columns)
            the_groups[the_groups.index.str.contains('GBM')] = 'GBM'

            keep = filter_cpm_by_group(
                the_data.divide(the_data.sum(), axis=1) * 1e6,
                the_groups,
                min_cpm=min_cpm
            )

            res_4.setdefault(m, {}).setdefault(pid, {})['iNSC'] = res_3[m][pid]['iNSC'].loc[keep]
            step4_n_filter.setdefault(m, {}).setdefault(pid, {})['iNSC'] = (keep.size, keep.sum())

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

    # report some of the numbers involved
    cols = []
    for k, nm in zip(['res_1', 'res_2', 'res_4'], ['original', 'original_filter', 'group_dispersion']):
        for m in methods:
            cols.append("%s_%s" % (nm, m))
    n_pair_only_intersect = pd.DataFrame(index=pids, columns=cols)

    for k, nm in zip(['res_1', 'res_2', 'res_4'], ['original', 'original_filter', 'group_dispersion']):
        for m in methods:
            this_res = to_save[k][m]

            # number of DE genes in Venn diagrams
            fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(10, 8))
            for i, pid in enumerate(pids):
                # number DE in the refs
                ax = axs.flat[i]
                venn.venn_diagram(*[this_res[pid][t].index for t in ['iNSC'] + refs], set_labels=['iNSC'] + refs, ax=ax)
                ax.set_title(pid, fontsize=16)
            for i in range(len(pids), 12):
                ax = axs.flat[i]
                ax.set_visible(False)
            fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.95)
            fig.savefig(os.path.join(outdir, "number_de_genes_ref_comparison_%s_%s.png" % (nm, m)), dpi=200)

            # number of PO genes in Venn diagrams
            fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(10, 6))
            for i, pid in enumerate(pids):
                a = this_res[pid]['iNSC'].index
                po = []
                for ref in refs:
                    b = this_res[pid][ref].index
                    vs, vc = setops.venn_from_arrays(a, b)
                    po.append(vs['10'])
                ax = axs.flat[i]
                venn.venn_diagram(*po, set_labels=refs, ax=ax)
                ax.set_title("GBM%s pair only" % pid, fontsize=16)
            for i in range(len(pids), 12):
                ax = axs.flat[i]
                ax.set_visible(False)
            fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.95)
            fig.savefig(os.path.join(outdir, "po_de_genes_ref_comparison_%s_%s.png" % (nm, m)), dpi=200)

            # number of PO DE genes

            # overlap between individual references in terms of PO genes shared
            pct_pair_only_intersect = pd.DataFrame(index=pids, columns=refs)

            for i, pid in enumerate(pids):
                a = this_res[pid]['iNSC'].index
                po = []
                for ref in refs:
                    b = this_res[pid][ref].index
                    vs, vc = setops.venn_from_arrays(a, b)
                    po.append(vs['10'])
                n_pair_only = np.array([len(t) for t in po])
                _, vc = setops.venn_from_arrays(*po)
                n_pair_only_intersect.loc[pid, "%s_%s" % (nm, m)] = vc['11']
                pct_pair_only_intersect.loc[pid] = float(vc['11']) / n_pair_only * 100.

            ax = pct_pair_only_intersect.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k', legend=False)
            ax.set_xlabel('Patient')
            ax.set_ylabel("% DE genes shared")
            ax.set_ylim([0, 100])
            ax.figure.tight_layout()
            ax.figure.savefig(os.path.join(outdir, "perc_po_gene_correspondence_%s_%s.png" % (nm, m)), dpi=200)

            plt.close('all')

    # effect of filtering on original approach
    cols_select = ['original_GLM', 'original_filter_GLM']
    legend_lbls = [
        'No filtering',
        'CPM filtering',
    ]
    to_plot = n_pair_only_intersect.loc[:, cols_select]
    to_plot.columns = legend_lbls
    ax = to_plot.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k')
    ax.set_ylabel("Number agreeing PO genes")
    ax.figure.savefig(os.path.join(outdir, "filtering_effect_number_agreeing_po.png"), dpi=200)

    # number of PO genes in the various different routines
    cols_select = ['original_filter_GLM', 'group_dispersion_GLM', 'group_dispersion_QLGLM']
    legend_lbls = [
        'Original, filter, GLM',
        'Group dispersion, GLM',
        'Group dispersion, QLGLM',
    ]
    to_plot = n_pair_only_intersect.loc[:, cols_select]
    to_plot.columns = legend_lbls
    ax = to_plot.plot.bar(color=['#ff9999', '#99cc99', '#9999ff'], ec='k')
    ax.set_ylabel("Number agreeing PO genes")
    ax.figure.savefig(os.path.join(outdir, "number_agreeing_po.png"), dpi=200)
