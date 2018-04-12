from rnaseq import loader, differential_expression, filter, general
from plotting import common, clustering
from stats import transformations
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
import seaborn as sns



def ecdf_func(x):
    x = np.sort(x)
    n = float(x.size)
    def result(v):
        return np.searchsorted(x, v, side='right') / n
    return result


def log_cpm(dat, base=2, offset=1.):
    dat = dat + offset
    if len(dat.shape) == 2:
        cpm = dat.divide(dat.sum(), axis=1) * 1e6
    else:
        cpm = dat.divide(dat.sum()) * 1e6
    return np.log(cpm) / np.log(base)


if __name__ == '__main__':
    pids = ['019', '031', '049', '052']
    min_cpm = 1

    ## 1) STAR CPM estimates

    ss2_obj = loader.load_references('wtchg_p180059', strandedness='u')
    assigned_sum = ss2_obj.data.sum()
    unassigned_sum = ss2_obj.data_unassigned.drop('N_unmapped').sum()

    ss2_pct_assigned = assigned_sum / (assigned_sum + unassigned_sum) * 100.

    print "SmartSeq2 samples % assigned"
    print ss2_pct_assigned

    polya_obj = loader.load_by_patient(pids)
    # restrict to relevant samples
    idx = (polya_obj.meta.type == 'iNSC')
    polya_obj.meta = polya_obj.meta.loc[idx]
    polya_obj.data = polya_obj.data.loc[:, polya_obj.meta.index]
    polya_obj.data_unassigned = polya_obj.data_unassigned.loc[:, polya_obj.meta.index]

    assigned_sum = polya_obj.data.sum()
    unassigned_sum = polya_obj.data_unassigned.drop('N_unmapped').sum()

    polya_pct_assigned = assigned_sum / (assigned_sum + unassigned_sum) * 100.

    print "Poly(A) samples % assigned"
    print polya_pct_assigned

    # combine data then eliminate genes that are not expressed
    ss2_data = ss2_obj.data
    ss2_data.columns = ["%s_smartseq" % t for t in ss2_data.columns]
    polya_data = polya_obj.data
    data = pd.concat((ss2_data, polya_data), axis=1)
    data = filter.filter_by_cpm(data, min_cpm=min_cpm, min_n_samples=1)

    # plot the CDF of the log2(CPM) values
    # for this purpose, we need to filter the CPM values for each col separately
    x_cdf = np.linspace(-5, 15, 500)
    log_cpm_ecdf = {}

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ss_nsc_lbl = False
    ss_opc_lbl = False
    pa_lbl = False

    for c, col in data.iteritems():
        this_dat = col.loc[col >= min_cpm] + 1
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = ecdf_func(this_log_cpm)
        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)
        log_cpm_ecdf[c] = this_log_cpm_ecdf

        lbl = None
        if 'smartseq' in c and 'OPC' in c:
            style = {'c': 'b'}
            if ss_opc_lbl is False:
                lbl = 'OPC SmartSeq2'
                ss_opc_lbl = True
        elif 'smartseq' in c and 'NSC' in c:
            style = {'c': 'r'}
            if ss_nsc_lbl is False:
                lbl = 'NSC SmartSeq2'
                ss_nsc_lbl = True
        else:
            style = {'c': 'k'}
            if pa_lbl is False:
                lbl = 'Poly(A)'
                pa_lbl = True
        ax.plot(x_cdf, this_log_cpm_ecdf, label=lbl, **style)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")

    # same, but only plot the SS2 samples, with labels included
    colours = common.COLOUR_BREWERS[ss2_obj.meta.shape[0]]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(ss2_data.columns):
        ax.plot(x_cdf, log_cpm_ecdf[c], c=colours[i], label=c)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")

    # same, but only plot the Poly(A) samples, with labels included
    colours = common.COLOUR_BREWERS[polya_obj.meta.shape[0]]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(polya_data.columns):
        ax.plot(x_cdf, log_cpm_ecdf[c], c=colours[i], label=c)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")

    # matching samples between two different preps

    matched_data = filter.filter_by_cpm(data.loc[:, data.columns.str.contains('NSC')], min_cpm=min_cpm, min_n_samples=1)
    matched_log_cpm = log_cpm(matched_data)
    row_colours= pd.DataFrame(common.COLOUR_BREWERS[2][0], index=matched_data.columns, columns=['Library'])
    row_colours.loc[row_colours.index.str.contains('smartseq')] = common.COLOUR_BREWERS[2][1]

    # clustering plot
    cg = clustering.plot_correlation_clustermap(matched_log_cpm, row_colors=row_colours)
    cg.gs.update(bottom=0.35, right=0.65)
    mad = transformations.median_absolute_deviation(matched_log_cpm).sort_values(ascending=False)
    cg = clustering.plot_correlation_clustermap(matched_log_cpm.loc[mad.index[:3000]], row_colors=row_colours)
    cg.gs.update(bottom=0.35, right=0.65)

    # each of the three pairings in each NSC
    for p in pids:
        this_matched = filter.filter_by_cpm(
            matched_data.loc[:, matched_data.columns.str.contains(p)],
            min_cpm=min_cpm,
            min_n_samples=1
        )
        this_matched_log_cpm = log_cpm(this_matched)
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
        ss_idx = this_matched.columns[this_matched.columns.str.contains('smartseq')]
        not_ss_idx = this_matched.columns[~this_matched.columns.str.contains('smartseq')]
        axs[0].scatter(this_matched_log_cpm.loc[:, not_ss_idx[0]], this_matched_log_cpm.loc[:, not_ss_idx[1]])
        axs[0].set_xlabel(not_ss_idx[0])
        axs[0].set_ylabel(not_ss_idx[1])
        axs[1].scatter(this_matched_log_cpm.loc[:, ss_idx[0]], this_matched_log_cpm.loc[:, not_ss_idx[0]])
        axs[1].set_xlabel(ss_idx[0])
        axs[1].set_ylabel(not_ss_idx[0])
        axs[2].scatter(this_matched_log_cpm.loc[:, ss_idx[0]], this_matched_log_cpm.loc[:, not_ss_idx[1]])
        axs[2].set_xlabel(ss_idx[0])
        axs[2].set_ylabel(not_ss_idx[1])

        fig.tight_layout()

    # run DE between the SmartSeq and regular samples
    de_method = 'QLGLM'

    # 1) ignoring patient, just consistent differences between SS and polyA
    patient = matched_data.columns.str.replace('_.*', '').str.replace('DURA', '')
    library_prep = ['SS' if 'smartseq' in t else 'polyA' for t in matched_data.columns]

    fit, design = differential_expression.edger_fit_glm(
        matched_data,
        de_method,
        '~patient + library_prep',
        patient=patient,
        library_prep=library_prep,
    )
    de_res_ss = differential_expression.edger_test(fit, design, "library_prepSS")
    general.add_gene_symbols_to_ensembl_data(de_res_ss)

    # 2) separate comparisons SS vs polyA
    de_res_separate = {}
    for p in pids:
        this_matched = filter.filter_by_cpm(
            matched_data.loc[:, matched_data.columns.str.contains(p)],
            min_cpm=min_cpm,
            min_n_samples=1
        )
        library_prep = ['SS' if 'smartseq' in t else 'polyA' for t in this_matched.columns]
        the_comparison = ('SS', 'polyA')
        de_res_separate[p] = differential_expression.run_one_de(this_matched, library_prep, the_comparison, method=de_method)


    # 3) how about DE between different NSC (within each library prep situation)?
    # this requires a comparison without replicates in the case of smartseq
    # must use GLM method here - QLGLM complains about estimating dof (?)
    de_method = 'GLM'
    comparisons_run = set()
    de_res_intra_ss = {}
    de_res_intra_pa = {}
    for p1 in pids:
        for p2 in pids:
            if p1 == p2:
                continue
            if (p1, p2) in comparisons_run:
                continue
            ss_idx = matched_data.columns.str.contains('smartseq') & matched_data.columns.str.contains("%s|%s" % (p1, p2))
            ss_data = matched_data.loc[:, ss_idx]
            ss_groups = patient[ss_idx]
            fit, design = differential_expression.edger_fit_glm(
                ss_data,
                de_method,
                '~0+patient',
                common_disp=True,
                patient=ss_groups,
            )
            de_res_intra_ss[(p1, p2)] = differential_expression.edger_test(fit, design, "patient%s - patient%s" % (p1, p2))

            # four ways to run the PA comparison (make it n=1 to be fair)
            pa_p1_idx = (~matched_data.columns.str.contains('smartseq')) & matched_data.columns.str.contains(p1)
            pa_p2_idx = (~matched_data.columns.str.contains('smartseq')) & matched_data.columns.str.contains(p2)
            pa_idx = pa_p1_idx | pa_p2_idx
            pa_data = matched_data.loc[:, pa_idx]
            pa_groups = patient[pa_idx]

            for i in range(0, 2):
                for j in range(0, 2):
                    this_p1 = matched_data.columns[pa_p1_idx][i]
                    this_p2 = matched_data.columns[pa_p2_idx][j]
                    this_pa_data = pa_data.loc[:, [this_p1, this_p2]]
                    this_pa_groups = pa_groups[pa_data.columns.isin([this_p1, this_p2])]

                    fit, design = differential_expression.edger_fit_glm(
                        this_pa_data,
                        de_method,
                        '~0+patient',
                        common_disp=True,
                        patient=this_pa_groups,
                    )
                    de_res_intra_pa[(this_p1, this_p2)] = differential_expression.edger_test(fit, design, "patient%s - patient%s" % (p1, p2))

            comparisons_run.add((p1, p2))
            comparisons_run.add((p2, p1))





