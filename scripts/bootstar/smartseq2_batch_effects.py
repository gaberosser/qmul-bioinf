from rnaseq import loader, differential_expression, filter, general
from plotting import common, clustering
from stats import transformations, basic
import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster import hierarchy as hc
import math
from matplotlib import pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from utils import output, setops
import references
import os


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
    min_cpm_individual = 0.1

    outdir = output.unique_output_dir("james_opc_smartseq2_vs_polya")

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

    # TMM normed version
    data_n = transformations.edger_tmm_normalisation_cpm(data)

    # plot the CDF of the log2(CPM) values
    # for this purpose, we need to filter the CPM values for each col separately
    x_cdf = np.linspace(-5, 15, 500)
    log_cpm_ecdf = {}

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ss_nsc_lbl = False
    ss_opc_lbl = False
    pa_lbl = False

    for c, this_dat in data.iteritems():
        this_dat = this_dat.loc[this_dat > 0]
        this_cpm = this_dat.divide(this_dat.sum()) * 1e6
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)
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
    fig.savefig(os.path.join(outdir, "ecdf_by_library_prep.png"), dpi=200)

    # again but TMM normed
    log_cpm_ecdf_tmm = {}


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ss_nsc_lbl = False
    ss_opc_lbl = False
    pa_lbl = False

    for c, this_cpm in data_n.iteritems():
        this_cpm = this_cpm.loc[this_cpm > 0]
        this_log_cpm = np.log2(this_cpm)
        this_ecdf_fun = basic.ecdf_func(this_log_cpm)
        this_log_cpm_ecdf = this_ecdf_fun(x_cdf)
        log_cpm_ecdf_tmm[c] = this_log_cpm_ecdf

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
    fig.savefig(os.path.join(outdir, "ecdf_by_library_prep_tmm.png"), dpi=200)

    # same, but only plot the SS2 samples, with labels included
    colours = common.COLOUR_BREWERS[ss2_obj.meta.shape[0]]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(ss2_data.columns):
        ax.plot(x_cdf, log_cpm_ecdf[c], c=colours[i], label=c)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")
    fig.savefig(os.path.join(outdir, "ecdf_smartseq2.png"), dpi=200)

    # TMM
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(ss2_data.columns):
        ax.plot(x_cdf, log_cpm_ecdf_tmm[c], c=colours[i], label=c)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")
    fig.savefig(os.path.join(outdir, "ecdf_smartseq2_tmm.png"), dpi=200)

    # same, but only plot the Poly(A) samples, with labels included
    colours = common.COLOUR_BREWERS[polya_obj.meta.shape[0]]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(polya_data.columns):
        ax.plot(x_cdf, log_cpm_ecdf[c], c=colours[i], label=c)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")
    fig.savefig(os.path.join(outdir, "ecdf_polya.png"), dpi=200)

    # TMM
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, c in enumerate(polya_data.columns):
        ax.plot(x_cdf, log_cpm_ecdf_tmm[c], c=colours[i], label=c)
    ax.legend(loc='lower right')
    ax.set_xlabel("log2(CPM)")
    ax.set_ylabel("Empirical CDF")
    fig.savefig(os.path.join(outdir, "ecdf_polya_tmm.png"), dpi=200)

    # matching samples between two different preps

    matched_data = filter.filter_by_cpm(data.loc[:, data.columns.str.contains('NSC')], min_cpm=min_cpm, min_n_samples=1)
    matched_log_cpm = log_cpm(matched_data)
    row_colours= pd.DataFrame(common.COLOUR_BREWERS[2][0], index=matched_data.columns, columns=['Library'])
    row_colours.loc[row_colours.index.str.contains('smartseq')] = common.COLOUR_BREWERS[2][1]

    # clustering plot
    cg = clustering.plot_correlation_clustermap(matched_log_cpm, row_colors=row_colours)
    cg.gs.update(bottom=0.35, right=0.65)
    cg.savefig(os.path.join(outdir, "cluster_log_cpm_corr_all_genes.png"), dpi=200)

    mad = transformations.median_absolute_deviation(matched_log_cpm).sort_values(ascending=False)
    cg = clustering.plot_correlation_clustermap(matched_log_cpm.loc[mad.index[:3000]], row_colors=row_colours)
    cg.gs.update(bottom=0.35, right=0.65)
    cg.savefig(os.path.join(outdir, "cluster_log_cpm_corr_3000_genes.png"), dpi=200)

    # repeat with TMM norming
    matched_log_cpm_n = transformations.edger_tmm_normalisation_cpm(matched_data)

    cg = clustering.plot_correlation_clustermap(matched_log_cpm_n, row_colors=row_colours)
    cg.gs.update(bottom=0.35, right=0.65)
    cg.savefig(os.path.join(outdir, "cluster_log_cpm_corr_all_genes_tmm.png"), dpi=200)

    mad = transformations.median_absolute_deviation(matched_log_cpm_n).sort_values(ascending=False)
    cg = clustering.plot_correlation_clustermap(matched_log_cpm_n.loc[mad.index[:3000]], row_colors=row_colours)
    cg.gs.update(bottom=0.35, right=0.65)
    cg.savefig(os.path.join(outdir, "cluster_log_cpm_corr_3000_genes_tmm.png"), dpi=200)

    # check with Spearman: in theory, TMM norming should make no difference?
    spearman_pdist = pd.DataFrame(index=matched_log_cpm.columns, columns=matched_log_cpm.columns, dtype=float)
    pairs_seen = set()
    for c1 in matched_log_cpm.columns:
        for c2 in matched_log_cpm.columns:
            if c1 == c2:
                spearman_pdist.loc[c1, c2] = 0
            if (c1, c2) in pairs_seen:
                continue
            pairs_seen.add((c1, c2))
            this_corr = stats.spearmanr(
                matched_log_cpm.loc[:, c1],
                matched_log_cpm.loc[:, c2]
            ).correlation
            spearman_pdist.loc[c1, c2] = 1 - this_corr
            spearman_pdist.loc[c2, c1] = 1 - this_corr
    X = spearman_pdist.values
    X[np.arange(len(X)), np.arange(len(X))] = 0
    spearman_cond_dist = hc.distance.squareform(X)

    cg = clustering.plot_correlation_clustermap(matched_log_cpm, row_colors=row_colours, distance=spearman_cond_dist)
    cg.gs.update(bottom=0.35, right=0.65)

    spearman_pdist_n = pd.DataFrame(index=matched_log_cpm_n.columns, columns=matched_log_cpm_n.columns, dtype=float)
    pairs_seen = set()
    for c1 in matched_log_cpm_n.columns:
        for c2 in matched_log_cpm_n.columns:
            if c1 == c2:
                spearman_pdist.loc[c1, c2] = 0
            if (c1, c2) in pairs_seen:
                continue
            pairs_seen.add((c1, c2))
            this_corr = stats.spearmanr(
                matched_log_cpm_n.loc[:, c1],
                matched_log_cpm_n.loc[:, c2]
            ).correlation
            spearman_pdist_n.loc[c1, c2] = 1 - this_corr
            spearman_pdist_n.loc[c2, c1] = 1 - this_corr
    X = spearman_pdist_n.values
    X[np.arange(len(X)), np.arange(len(X))] = 0
    spearman_cond_dist_n = hc.distance.squareform(X)

    cg = clustering.plot_correlation_clustermap(matched_log_cpm_n, row_colors=row_colours, distance=spearman_cond_dist_n)
    cg.gs.update(bottom=0.35, right=0.65)

    # each of the three pairings in each NSC
    for p in pids:
        this_matched = filter.filter_by_cpm(
            matched_data.loc[:, matched_data.columns.str.contains(p)],
            min_cpm=min_cpm,
            min_n_samples=1
        )
        this_matched_log_cpm = log_cpm(this_matched)
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(11, 3.5))
        ss_idx = this_matched.columns[this_matched.columns.str.contains('smartseq')]
        not_ss_idx = this_matched.columns[~this_matched.columns.str.contains('smartseq')]

        a = this_matched_log_cpm.loc[:, not_ss_idx[0]]
        b = this_matched_log_cpm.loc[:, not_ss_idx[1]]
        axs[0].scatter(a, b)
        axs[0].set_xlabel(not_ss_idx[0])
        axs[0].set_ylabel(not_ss_idx[1])
        axs[0].text(-5, 12.5, "Spearman rho: %.3f" % stats.spearmanr(a, b).correlation)

        a = this_matched_log_cpm.loc[:, ss_idx[0]]
        b = this_matched_log_cpm.loc[:, not_ss_idx[0]]
        axs[1].scatter(a, b)
        axs[1].set_xlabel(ss_idx[0])
        axs[1].set_ylabel(not_ss_idx[0])
        axs[1].text(-5, 12.5, "Spearman rho: %.3f" % stats.spearmanr(a, b).correlation)

        a = this_matched_log_cpm.loc[:, ss_idx[0]]
        b = this_matched_log_cpm.loc[:, not_ss_idx[1]]
        axs[2].scatter(a, b)
        axs[2].set_xlabel(ss_idx[0])
        axs[2].set_ylabel(not_ss_idx[1])
        axs[2].text(-5, 12.5, "Spearman rho: %.3f" % stats.spearmanr(a, b).correlation)

        for i in range(3):
            axs[i].set_xlim([-6, 15])
            axs[i].set_ylim([-6, 15])
            axs[i].set_aspect('equal')

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "correlation_plots_%s.png" % p), dpi=200)

    # run DE between the SmartSeq and regular samples
    de_method = 'QLGLM'

    # 1) ignoring patient, just consistent differences between SS and polyA
    # in light of the 3 subgroups in the ECDF, this probably contains too many different effects to be useful
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
    # restrict this to matching passage / subclone
    # use grouped dispersion approach (SS vs PolyA) and GLM method to avoid error when estimating DOF
    de_method = 'GLM'

    sample_pairs = {
        '019': ('DURA019_NSC_N8C_P2', 'DURA019_NSC_N8C_smartseq'),
        '031': ('DURA031_NSC_N44B_P2', 'DURA031_NSC_N44B_smartseq'),
        '049': ('DURA049_NSC_N19_P4', 'DURA049_NSC_N19_P4_smartseq'),
        '052': ('DURA052_NSC_N4_P3', 'DURA052_NSC_N4_P3_smartseq')
    }
    all_snames = []
    [all_snames.extend(sample_pairs[p]) for p in pids]
    this_matched = matched_data.loc[:, all_snames]

    groups = pd.Series(index=this_matched.columns)
    for p in pids:
        sp = sample_pairs[p]
        groups.loc[sp[0]] = "PolyA%s" % p
        groups.loc[sp[1]] = "SmartSeq%s" % p

    groups_for_disp = pd.Series('PolyA', index=this_matched.columns)
    groups_for_disp.loc[groups_for_disp.index.str.contains('smartseq')] = 'SmartSeq'

    # 2a) Pool across library prep methods for dispersion estimation
    dgelist, _ = differential_expression.edger_dgelist(
        this_matched,
        the_formula="~0 + groups",
        groups=groups_for_disp.values
    )
    design = differential_expression.model_matrix_from_formula("~0 + group", group=groups.values)
    fit, _ = differential_expression.edger_fit_dgelist(dgelist, de_method, design=design)

    de_res_separate = {}
    for p in pids:
        de_res_separate[p] = differential_expression.edger_test(fit, design, "groupSmartSeq%s - groupPolyA%s" % (p, p))
        general.add_gene_symbols_to_ensembl_data(de_res_separate[p])
        print "Patient %s: %d DE genes in SmartSeq2 - PolyA (%d up, %d down)." % (
            p,
            de_res_separate[p].shape[0],
            (de_res_separate[p].logFC > 0).sum(),
            (de_res_separate[p].logFC < 0).sum(),
        )

    de_in_all = references.ensembl_to_gene_symbol(
        setops.reduce_intersection(*[t.index for t in de_res_separate.values()])
    )

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10, 10))
    fig.tight_layout()
    for i, p in enumerate(pids):
        ax = axs.flat[i]
        sp = sample_pairs[p]

        for_scatter = []
        for s in sp:
            this_dat = this_matched.loc[de_res_separate[p].index, s] + 1
            this_cpm = this_dat.divide(this_dat.sum()) * 1e6
            this_log_cpm = np.log2(this_cpm)
            for_scatter.append(this_log_cpm)

        abs_fc = np.abs(de_res_separate[p].logFC).sort_values(ascending=False)

        ix_up = de_res_separate[p].index[de_res_separate[p].logFC > 0]
        ix_down = de_res_separate[p].index[de_res_separate[p].logFC < 0]

        ax.scatter(
            for_scatter[0].loc[ix_up],
            for_scatter[1].loc[ix_up],
            color='#77d86c'
        )
        ax.scatter(
            for_scatter[0].loc[ix_down],
            for_scatter[1].loc[ix_down],
            color='#9ed3ff'
        )
        ax.set_aspect('equal')
        ax.plot([-1, 20], [-1, 20], 'k--')

        # annotate the largest FC genes (top 10%)
        to_annot = abs_fc.index[:int(np.ceil(0.1 * abs_fc.shape[0]))]
        to_annot = references.ensembl_to_gene_symbol(to_annot)
        # to_annot.loc[to_annot.isnull()] = to_annot.index[to_annot.isnull()]
        to_annot = to_annot.dropna()
        texts = []
        for e, txt in to_annot.iteritems():
            weight = None
            if e in de_in_all:
                weight = 'bold'
            if e in ix_up:
                c = 'g'
            else:
                c = 'b'

            tt = ax.text(
                for_scatter[0].loc[e],
                for_scatter[1].loc[e],
                txt,
                color=c,
                weight=weight
            )
            texts.append(tt)
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'), ax=ax)
    fig.savefig(os.path.join(outdir, "scatter_logcpm_de_genes.png"), dpi=200)

    for i, p in enumerate(pids):
        sp = sample_pairs[p]
        fig = plt.figure()
        ax = fig.add_subplot(111)

        cs = ['b', 'g']
        ls = ['%s PolyA' % p, '%s SmartSeq2' % p]
        ecdf_funs = []
        log_cpms = []

        for j, s in enumerate(sp):
            this_dat = this_matched[s].loc[this_matched[s] >= min_cpm] + 1
            this_cpm = this_dat.divide(this_dat.sum()) * 1e6
            this_log_cpm = np.log2(this_cpm)
            log_cpms.append(this_log_cpm)
            this_ecdf_fun = basic.ecdf_func(this_log_cpm)
            ecdf_funs.append(this_ecdf_fun)
            this_log_cpm_ecdf = this_ecdf_fun(x_cdf)
            ax.plot(x_cdf, this_log_cpm_ecdf, c=cs[j], label=ls[j])

        # now find the ECDF values at the logCPM values of genes of interest and plot
        abs_fc = np.abs(de_res_separate[p].logFC).sort_values(ascending=False)
        fc = de_res_separate[p].logFC.loc[abs_fc.index]

        for e in abs_fc.index:
            try:
                lc0 = log_cpms[0].loc[e]
                ec0 = ecdf_funs[0](lc0)
            except KeyError:
                lc0 = None
                ec0 = None
            try:
                lc1 = log_cpms[1].loc[e]
                ec1 = ecdf_funs[1](lc1)
            except KeyError:
                lc1 = None
                ec1 = None
            if lc0 is None and lc1 is None:
                print "Something is wrong here: can't find gene %s in either of the arrays for PID %s" % (e, p)
            elif lc0 is not None and lc1 is not None:
                if lc0 > lc1:
                    # higher expression in poly(A): blue line, linked upwards
                    this_max_y = max(ec0, ec1)
                    ax.plot([lc1, lc1], [ec1, this_max_y], c='b', alpha=0.2)
                    ax.plot([lc1, lc0], [this_max_y, this_max_y], c='b', alpha=0.2)
                else:
                    # higher expression in SS2: green line, linked downwards
                    this_min_y = min(ec0, ec1)
                    ax.plot([lc1, lc1], [ec1, this_min_y], c='g', alpha=0.2)
                    ax.plot([lc1, lc0], [this_min_y, this_min_y], c='g', alpha=0.2)
        ax.set_title("%s poly(A) vs SmartSeq2" % p)
        ax.set_xlabel('log2(CPM)')
        ax.set_ylabel('ECDF')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'ecdf_with_de_interactions_%s.png' % p), dpi=200)

    # 2b) Pool across all samples for dispersion estimation
    # this common dispersion approach is WAY too conservative as it aggregates ALL samples for dispersion estimation
    # (resulting in a very high estimate)
    fit, design = differential_expression.edger_fit_glm(
        this_matched,
        de_method,
        '~0+group',
        group=groups.values,
        common_disp=True
    )
    de_res_separate_common = {}
    for p in pids:
        de_res_separate_common[p] = differential_expression.edger_test(fit, design, "groupSmartSeq%s - groupPolyA%s" % (p, p))
        general.add_gene_symbols_to_ensembl_data(de_res_separate_common[p])


    # 3) how about DE between different NSC (within each library prep situation)?
    # this requires a comparison without replicates in the case of smartseq
    # must use GLM method here - QLGLM complains about estimating dof
    de_method = 'GLM'

    # 3a) Matching the format of 2a: 

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





