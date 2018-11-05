import os
from settings import HGIC_LOCAL_DIR
import pandas as pd
from utils import output, setops, excel, ipa
from plotting import common
import collections
import consts
from matplotlib import pyplot as plt, patches, collections, gridspec
import seaborn as sns
import numpy as np
from scipy import stats
from stats import nht, basic


def ipa_results_to_wideform(res, plogalpha):
    """
    Convert the IPA results dictionary into a wideform pd.DataFrame.
    Owing to the potentially large number of comparisons, we can't use the Venn approach here, but there's no need.
    :param res:
    :param plogalpha:
    :return:
    """
    de_all_pathways = sorted(setops.reduce_union(*[t.index for t in res.values()]))
    export_wideform = pd.DataFrame(index=de_all_pathways)
    member_cols = []
    for k, v in res.items():
        sign_ix = v.index[v['-logp'] >= plogalpha]
        this_yn = pd.Series('N', index=de_all_pathways)
        this_yn.loc[sign_ix] = 'Y'
        member_cols.append(k)
        export_wideform.insert(
            export_wideform.shape[1],
            k,
            this_yn
        )
        for col in ['-logp', 'z', 'ratio', 'n_gene']:
            export_wideform.insert(
                export_wideform.shape[1],
                "%s_%s" % (k, col),
                v.reindex(de_all_pathways)[col]
            )

    # add n gene in pathway as single const column
    rr = export_wideform.loc[:, export_wideform.columns.str.contains('ratio')]
    ng = export_wideform.loc[:, export_wideform.columns.str.contains('n_gene')]
    n_gene_tot = (ng.astype(float).values / rr.astype(float)).mean(axis=1).round().astype(int)
    export_wideform.insert(0, 'n_gene_in_pathway', n_gene_tot)

    return export_wideform


def generate_summary_df(all_in, comps, pids=consts.PIDS):
    """
    Given the full results dictionary, generate a single DataFrame summary of the number of comparisons and IDs
    exhibiting each pathway.
    :param all_in: DataFrame, indexes are pathways and columns are comparison names. Entries are boolean, where a T
    indicates that the pathway is significant in that comparison
    :param comps: List of comparisons, not including syngeneic
    :param pids: List of PIDs.
    :return:
    """
    pathways = all_in.index

    n_set = pd.DataFrame(
        0,
        index=pathways,
        columns=['Syngen. only', 'Ref. only', 'Intersect.'],
        dtype=int
    )
    so = dict([(pw, []) for pw in pathways])
    ro = dict([(pw, []) for pw in pathways])
    inters = dict([(pw, []) for pw in pathways])
    for pid in pids:
        s = all_in.index[all_in.loc[:, "%s_syngeneic" % pid]]
        r = all_in.index[all_in.loc[:, ["%s_%s" % (pid, t) for t in comps]].any(axis=1)]
        vs, _ = setops.venn_from_arrays(s, r)

        n_set.loc[vs['10'], 'Syngen. only'] += 1
        n_set.loc[vs['01'], 'Ref. only'] += 1
        n_set.loc[vs['11'], 'Intersect.'] += 1

        for pw in vs['10']:
            so[pw].append(pid)
        for pw in vs['01']:
            ro[pw].append(pid)
        for pw in vs['11']:
            inters[pw].append(pid)

    # output excel file giving at-a-glance access to which patients are involved in each pathway, categorised as
    # 'syn only', 'ref only' and 'intersection'
    at_a_glance = pd.DataFrame(
        index=pathways,
        columns=['n_syngen_only', 'syngen_only_pids', 'n_ref_only', 'ref_only_pids', 'n_intersect', 'intersect_pids'],
        dtype=object
    )
    for pw in pathways:
        at_a_glance.loc[pw, 'n_syngen_only'] = len(so[pw])
        at_a_glance.loc[pw, 'syngen_only_pids'] = ';'.join(so[pw])
        at_a_glance.loc[pw, 'n_ref_only'] = len(ro[pw])
        at_a_glance.loc[pw, 'ref_only_pids'] = ';'.join(ro[pw])
        at_a_glance.loc[pw, 'n_intersect'] = len(inters[pw])
        at_a_glance.loc[pw, 'intersect_pids'] = ';'.join(inters[pw])

    return n_set, at_a_glance


def compute_enrichment_combined_references(df):
    """
    Compute the effect and significance of enrichment of pathways in either syngeneic or reference comparisons.
    :param df: DataFrame containing two columns, 'syn' and 'ref'. These can either be the count (number of
    comparisons identifying this pathway) or the sum of -logp values ('total significance' of this pathway).
    :return:
    """
    N = df.shape[0]
    # Ntot is the number of pairwise comparisons
    Ntot = N * (N - 1) / 2.

    this_wsrt = nht.wilcoxon_signed_rank_statistic(df.syn, df.ref, zero_method='pratt')
    this_p = nht.wilcoxon_signed_rank_test(df.syn, df.ref, distribution='exact')
    this_rank_corr = (this_wsrt['r_plus'] - this_wsrt['r_minus']) / Ntot

    print "Comparing the NUMBER of samples showing enrichment in a given pathway, combining references."
    if this_p < 0.05:
        print "Reject null (p=%.3f). Effect size/direction: %.3f (%s)." % (
            this_p,
            this_rank_corr,
            "syn > ref" if this_rank_corr > 0 else "ref > syn"
        )

    res = {
        'pval': this_p,
        'rank_correlation': this_rank_corr,
    }
    res.update(this_wsrt)
    return res


def compute_enrichment_separate_references(df, comps):
    """
    Compute the relative enrichment of pathways identified by the syngeneic comparison, vs each of the other
    comparisons.
    :param df: DataFrame, index is pathway, column is comparison. Value is whatever we want to test (e.g. count,
    sum plogp).
    :param comps: Array of comparison names. The first should be 'syngeneic' (or whatever we are testing).
    :return:
    """
    N = df.shape[0]
    # Ntot is the number of pairwise comparisons
    Ntot = N * (N - 1) / 2.

    wsrt = {}
    pval = {}
    rank_corr = {}

    c_base = comps[0]

    for c in comps[1:]:
        wsrt[c] = nht.wilcoxon_signed_rank_statistic(df[c_base], df[c], zero_method='pratt')
        pval[c] = nht.wilcoxon_signed_rank_test(df[c_base], df[c], distribution='exact')
        rank_corr[c] = (wsrt[c]['r_plus'] - wsrt[c]['r_minus']) / Ntot

        if pval[c] < 0.05:
            print "Comparison syngeneic vs %s. Reject null (p=%.3f). Effect size/direction: %.3f (%s)." % (
                c,
                pval[c],
                rank_corr[c],
                "syn > ref" if rank_corr[c] > 0 else "ref > syn"
            )
        else:
            print "Comparison syngeneic vs %s. Do not reject null (p=%.3f)." % (
                c,
                pval[c]
            )

    res = {
        'p': pval,
        'rank_correlation': rank_corr,
        'statistic': wsrt
    }
    return res


def generate_plotting_structures(res, pathways, plogalpha):
    all_p = {}
    all_z = {}
    all_in = {}

    for k, v in res.items():
        # only keep significant pathways
        the_dat = v.reindex(pathways)

        all_in[k] = the_dat['-logp'] > plogalpha
        all_p[k] = the_dat['-logp'].mask(the_dat['-logp'] < plogalpha)
        all_z[k] = the_dat['z']

    all_in = pd.DataFrame(all_in).fillna(False)
    all_z = pd.DataFrame(all_z)
    all_p = pd.DataFrame(all_p)

    return {
        'in': all_in,
        'z': all_z,
        'p': all_p
    }


def pathway_involvement_heatmap_by_p(
        p_dat,
        n_set,
        pathway_order,
        pids,
        comparison_dict,
        orientation='vertical',
        vmin=None,
        vmax=None,
        count_cmap='Blues'
):
    n = len(pids) + 2
    m = 5

    if orientation == 'vertical':
        figsize = (9., 11.5)
        gs_kw = dict(
            left=0.43,
            right=0.93,
            top=0.998,
            bottom=0.1,
            wspace=0.1,
            hspace=0.1,
        )
        ncols = n
        nrows = m
    else:
        figsize = (14, 9.)
        ncols = m
        nrows = n
        raise NotImplementedError("TODO: define gs_kw.")

    fig = plt.figure(figsize=figsize)

    # set up axis grid
    gs = gridspec.GridSpec(
        ncols=ncols,
        nrows=nrows,
        width_ratios=[1] * len(pids) + [1, .5],
        **gs_kw
    )
    if orientation == 'vertical':
        cax = fig.add_subplot(gs[1:-1, -1])
    else:
        cax = fig.add_subplot(gs[-1, 1:-1])

    axs = []

    for i, pid in enumerate(pids):
        if orientation == 'vertical':
            ax = fig.add_subplot(gs[:, i])
        else:
            ax = fig.add_subplot(gs[i, :])

        axs.append(ax)
        this_dat = []
        for j, c in enumerate(sorted(comparison_dict.keys())):
            # embed this set of data in full sorted list
            the_dat = p_dat.loc[:, '%s_%s' % (pid, c)].reindex(pathway_order)
            this_dat.append(the_dat)

        if orientation == 'vertical':
            this_dat = pd.concat(this_dat, axis=1).iloc[:, ::-1]
        else:
            ## TODO: check this
            this_dat = pd.concat(this_dat, axis=1).transpose()

        sns.heatmap(
            this_dat,
            mask=this_dat.isnull(),
            cmap='YlOrRd',
            linewidths=.2,
            linecolor='w',
            vmin=vmin,
            vmax=vmax,
            yticklabels=False,
            cbar=i == 0,
            cbar_ax=cax,
            cbar_kws={"orientation": orientation},
            ax=ax,
        )
        if orientation == 'vertical':
            cax.set_ylabel('$-\log(p)$')
            ax.xaxis.set_ticklabels(comparison_dict.values(), rotation=90.)
            ax.set_xlabel(pid)
        else:
            cax.set_xlabel('$-\log(p)$')
            ax.yaxis.set_ticklabels(comparison_dict.values(), rotation=0.)
            ax.set_ylabel(pid)

        ax.set_facecolor('0.6')

    eax = fig.add_subplot(gs[:, len(pids)])
    # need to leave ticklabels on, or annotations don't show?
    if orientation == 'vertical':
        n_set_plot = n_set.loc[pathway_order]
    else:
        n_set_plot = n_set.loc[pathway_order].transpose()
    sns.heatmap(
        n_set_plot,
        mask=n_set_plot == 0,
        cmap=count_cmap,
        cbar=False,
        ax=eax,
        annot=True,
        fmt="d",
        annot_kws={"size": 6}
    )
    if orientation == 'vertical':
        eax.yaxis.set_ticks([])
        plt.setp(eax.xaxis.get_ticklabels(), rotation=90)
        axs[0].set_yticks(0.5 + np.arange(len(pathway_order)))
        axs[0].set_yticklabels(pathway_order[::-1], rotation=0, fontsize=7)
    else:
        eax.xaxis.set_ticks([])
        plt.setp(eax.yaxis.get_ticklabels(), rotation=0)
        # TODO: check this is correct
        axs[-1].set_xticks(0.5 + np.arange(len(pathway_order)))
        axs[-1].set_xticklabels(pathway_order, rotation=90, fontsize=7)

    # colorbar outline
    cbar = axs[0].collections[0].colorbar
    cbar.outline.set_linewidth(1.)
    cbar.outline.set_edgecolor('k')

    return {
        'figure': fig,
        'axs': axs,
        'cax': cax,
        'cbar': cbar,
        'eax': eax,
    }


def plot_delta_histogram(
    y,
    x,
    nbin=10,
    fmin=0.01,
    fmax=0.99,
    ax=None,
    summary_metric='median'
):
    delta = y - x
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    dsort = sorted(delta)
    dmin = dsort[int(fmin * len(dsort))]
    dmax = dsort[int(fmax * len(dsort))]

    edges = np.linspace(dmin, dmax, nbin)
    # make this symmetric about zero
    if (dmin < 0) and (dmax > 0):
        ix = np.searchsorted(edges, 0, side='left')
        dx = edges[1] - edges[0]
        edges += -edges[ix] - dx / 2.

    ax.hist(delta, edges)
    ax.set_xlabel('Difference in number of comparisons (syn - ref)')
    ax.axvline(getattr(delta, summary_metric)(), ls='--', color='k', label=summary_metric)

    return ax


if __name__ == '__main__':
    # set a minimum pval for pathways to be used
    alpha = 0.005
    plogalpha = -np.log10(alpha)
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    de_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways')
    dm_indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/methylation/merged_s1_s2/ipa/pathways')
    outdir = output.unique_output_dir()

    # keys are the term used in the filename, values are those used in the columns
    ref_names = ['h9', 'gibco']
    comps = {
        'syngeneic': 'syngeneic',
        'h9': 'H9',
        'gibco': 'GIBCO'
    }
    comparison_names = {
        'syngeneic': 'Syngen.',
        'h9': 'H9',
        'gibco': 'Gibco'
    }
    pids = consts.PIDS

    #######################################################
    # DE
    #######################################################

    # load IPA results from raw data and combine into a single export file
    # format: Excel, wideform
    de_res = ipa.load_raw_reports(
        de_indir,
        'de_s2_{0}_{1}.txt',
        pids,
        comps
    )
    for k, v in de_res.items():
        rele_ix = v.index[v['-logp'] >= plogalpha_relevant]
        de_res['_'.join(k)] = v.loc[rele_ix]
        de_res.pop(k)

    # wideform version of this (i.e. 30 blocks)
    de_res_wide = ipa_results_to_wideform(de_res, plogalpha)
    de_res_wide.to_excel(os.path.join(outdir, "full_de_ipa_results.xlsx"))

    # get a list of significant pathways (in at least one comparison)
    de_pathways_significant = set()
    for k, v in de_res.items():
        de_pathways_significant.update(v.index[v['-logp'] > plogalpha])

    # export significant results to an Excel file with separate tabs
    de_res_sign = dict([
        (k, v.loc[v['-logp'] > plogalpha]) for k, v in de_res.items()
    ])
    excel.pandas_to_excel(de_res_sign, os.path.join(outdir, "full_de_ipa_results_significant_separated.xlsx"))

    # export wideform, reduced to include only significant pathways
    de_res_wide.loc[sorted(de_pathways_significant)].to_excel(
        os.path.join(outdir, "full_de_ipa_results_significant.xlsx")
    )

    # useful structures for plots
    dd = generate_plotting_structures(de_res, de_pathways_significant, plogalpha)
    de_all_p = dd['p']
    de_all_z = dd['z']
    de_all_in = dd['in']
    
    # at-a-glance export
    de_n_set, at_a_glance = generate_summary_df(de_all_in, ref_names)
    at_a_glance.to_excel(os.path.join(outdir, "de_ipa_results_patients_by_s2_category.xlsx"))

    # plot 1) P values, ordered by sum of -log(P)
    p_order = de_all_p.sum(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        de_all_p,
        de_n_set,
        p_order,
        pids,
        comparison_names
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.tiff"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.pdf"), dpi=200)

    # plot 2) P values, ordered by mean of -log(P)
    p_order = de_all_p.mean(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        de_all_p,
        de_n_set,
        p_order,
        pids,
        comparison_names
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.tiff"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.pdf"), dpi=200)

    # Test relative enrichment of pathway detection by syngeneic and references
    # 1a) Combined references, number of detections

    nn = de_n_set.iloc[:, :2].add(de_n_set.iloc[:, -1], axis=0)
    nn.columns = ['syn', 'ref']
    print "DE. COMBINED references. Number showing enrichment in a given pathway."
    wsrt_res = compute_enrichment_combined_references(nn)

    ax = plot_delta_histogram(nn.syn, nn.ref, nbin=10)
    ax.figure.savefig(os.path.join(outdir, "de_histogram_delta_number_comparisons_references_together.png"), dpi=200)

    # 1b) Combined references, sum of plogp

    pp = pd.DataFrame(index=de_pathways_significant, columns=['syn', 'ref'])
    pp.loc[:, 'syn'] = de_all_p[de_all_p.columns[de_all_p.columns.str.contains('syngeneic')]].sum(axis=1)
    pp.loc[:, 'ref'] = de_all_p[de_all_p.columns[~de_all_p.columns.str.contains('syngeneic')]].sum(axis=1)
    print "DE. COMBINED references. Sum of plogp showing enrichment in a given pathway."
    wsrt_res = compute_enrichment_combined_references(pp)

    ax = plot_delta_histogram(pp.syn, pp.ref, nbin=20)
    ax.figure.savefig(os.path.join(outdir, "de_histogram_delta_sum_plogp_comparisons_references_together.png"), dpi=200)

    # 2a) Consider references separately, number of detections
    nn = {}
    for c in comps:
        this = de_all_in.loc[:, de_all_in.columns.str.contains(c)]
        nn[c] = this.sum(axis=1)
    nn = pd.DataFrame(nn)
    print "DE. SEPARATE references. Number showing enrichment in a given pathway"
    wsrt_res = compute_enrichment_separate_references(nn, ['syngeneic'] + ref_names)

    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)
    i = 0
    for c in ref_names:
        plot_delta_histogram(nn['syngeneic'], nn[c], nbin=10, ax=axs[i], summary_metric='mean')
        axs[i].set_xlabel("# syngeneic - # %s" % c.title())
        i += 1
    axs[-1].legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_histogram_delta_number_comparisons.png"), dpi=200)

    # plot ECDF (ish)
    p_order = de_all_in.sum(axis=1).sort_values(ascending=False).index
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for c in comps:
        ax.plot(nn[c].loc[p_order].values.cumsum(), label=c.title())
    ax.set_xlabel('Ranked pathway')
    ax.set_ylabel('Cumulative number of patients with enrichment')
    ax.legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_number_comparisons_ranked_cumul_sum.png"), dpi=200)

    # 2b) Consider references separately, sum of -logp
    pp = {}
    for c in comps:
        this = de_all_p.loc[:, de_all_p.columns.str.contains(c)]
        pp[c] = this.sum(axis=1)
    pp = pd.DataFrame(pp)
    print "DE. SEPARATE references. Sum of -logp for each given pathway"
    wsrt_res = compute_enrichment_separate_references(pp, ['syngeneic'] + ref_names)

    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)
    i = 0
    for c in ref_names:
        plot_delta_histogram(pp['syngeneic'], pp[c], nbin=20, ax=axs[i])
        axs[i].set_xlabel("# syngeneic - # %s" % c.title())
        i += 1

    axs[-1].legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_histogram_delta_sum_pvalues.png"), dpi=200)

    # plot ECDF (ish)
    p_order = de_all_p.sum(axis=1).sort_values(ascending=False).index
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for c in comps:
        ax.plot(pp[c].loc[p_order].values.cumsum(), label=c.title())
    ax.set_xlabel('Ranked pathway')
    ax.set_ylabel('Cumulative sum of -log10(p)')
    ax.legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "de_sum_pvalues_comparisons_ranked_cumul_sum.png"), dpi=200)

    #######################################################
    # DMR
    #######################################################

    # load IPA results from raw data and combine into a single export file
    # format: Excel, wideform
    dm_res = ipa.load_raw_reports(
        dm_indir,
        'dmr_s2_{0}_{1}.txt',
        pids,
        comps
    )
    for k, v in dm_res.items():
        rele_ix = v.index[v['-logp'] >= plogalpha_relevant]
        dm_res['_'.join(k)] = v.loc[rele_ix]
        dm_res.pop(k)

    # wideform version of this (i.e. 30 blocks)
    dm_res_wide = ipa_results_to_wideform(dm_res, plogalpha)
    dm_res_wide.to_excel(os.path.join(outdir, "full_dmr_ipa_results.xlsx"))

    # get a list of significant pathways (in at least one comparison)
    dm_pathways_significant = set()
    for k, v in dm_res.items():
        dm_pathways_significant.update(v.index[v['-logp'] > plogalpha])

    # export significant results to an Excel file with separate tabs
    dm_res_sign = dict([
        (k, v.loc[v['-logp'] > plogalpha]) for k, v in dm_res.items()
    ])
    excel.pandas_to_excel(dm_res_sign, os.path.join(outdir, "full_dmr_ipa_results_significant_separated.xlsx"))

    # export wideform, reduced to include only significant pathways
    dm_res_wide.loc[sorted(dm_pathways_significant)].to_excel(
        os.path.join(outdir, "full_dmr_ipa_results_significant.xlsx")
    )

    # useful structures for plots
    dd = generate_plotting_structures(dm_res, dm_pathways_significant, plogalpha)
    dm_all_p = dd['p']
    dm_all_z = dd['z']
    dm_all_in = dd['in']

    # at-a-glance export
    dm_n_set, at_a_glance = generate_summary_df(dm_all_in, ref_names)
    at_a_glance.to_excel(os.path.join(outdir, "dmr_ipa_results_patients_by_s2_category.xlsx"))

    # plot 1) P values, ordered by sum of -log(P)
    p_order = dm_all_p.sum(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        dm_all_p,
        dm_n_set,
        p_order,
        pids,
        comparison_names,
        count_cmap='Greens',
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_dmr.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_dmr.tiff"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_dmr.pdf"), dpi=200)

    # plot 2) P values, ordered by sum of -log(P)
    p_order = dm_all_p.mean(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        dm_all_p,
        dm_n_set,
        p_order,
        pids,
        comparison_names,
        count_cmap='Greens',
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_dmr.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_dmr.tiff"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_dmr.pdf"), dpi=200)

    # Test relative enrichment of pathway detection by syngeneic and references
    # 1a) Combined references, number of detections

    nn = dm_n_set.iloc[:, :2].add(dm_n_set.iloc[:, -1], axis=0)
    nn.columns = ['syn', 'ref']
    print "DMR. COMBINED references. Number showing enrichment in a given pathway."
    wsrt_res = compute_enrichment_combined_references(nn)

    ax = plot_delta_histogram(nn.syn, nn.ref, nbin=10)
    ax.figure.savefig(os.path.join(outdir, "dmr_histogram_delta_number_comparisons_references_together.png"), dpi=200)

    # 1b) Combined references, sum of plogp

    pp = pd.DataFrame(index=dm_pathways_significant, columns=['syn', 'ref'])
    pp.loc[:, 'syn'] = dm_all_p[dm_all_p.columns[dm_all_p.columns.str.contains('syngeneic')]].sum(axis=1)
    pp.loc[:, 'ref'] = dm_all_p[dm_all_p.columns[~dm_all_p.columns.str.contains('syngeneic')]].sum(axis=1)
    print "DMR. COMBINED references. Sum of plogp showing enrichment in a given pathway."
    wsrt_res = compute_enrichment_combined_references(pp)

    ax = plot_delta_histogram(pp.syn, pp.ref, nbin=20)
    ax.figure.savefig(os.path.join(outdir, "dmr_histogram_delta_sum_plogp_comparisons_references_together.png"), dpi=200)

    # 2a) Consider references separately, number of detections
    nn = {}
    for c in comps:
        this = dm_all_in.loc[:, dm_all_in.columns.str.contains(c)]
        nn[c] = this.sum(axis=1)
    nn = pd.DataFrame(nn)
    print "DMR. SEPARATE references. Number showing enrichment in a given pathway"
    wsrt_res = compute_enrichment_separate_references(nn, ['syngeneic'] + ref_names)

    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)
    i = 0
    for c in ref_names:
        plot_delta_histogram(nn['syngeneic'], nn[c], nbin=10, ax=axs[i], summary_metric='mean')
        axs[i].set_xlabel("# syngeneic - # %s" % c.title())
        i += 1
    axs[-1].legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_histogram_delta_number_comparisons.png"), dpi=200)

    # plot ECDF (ish)
    p_order = dm_all_in.sum(axis=1).sort_values(ascending=False).index
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for c in comps:
        ax.plot(nn[c].loc[p_order].values.cumsum(), label=c.title())
    ax.set_xlabel('Ranked pathway')
    ax.set_ylabel('Cumulative number of patients with enrichment')
    ax.legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_number_comparisons_ranked_cumul_sum.png"), dpi=200)

    # 2b) Consider references separately, sum of -logp
    pp = {}
    for c in comps:
        this = dm_all_p.loc[:, dm_all_p.columns.str.contains(c)]
        pp[c] = this.sum(axis=1)
    pp = pd.DataFrame(pp)
    print "DMR. SEPARATE references. Sum of -logp for each given pathway"
    wsrt_res = compute_enrichment_separate_references(pp, ['syngeneic'] + ref_names)

    fig, axs = plt.subplots(ncols=2, sharex=True, sharey=True)
    i = 0
    for c in ref_names:
        plot_delta_histogram(pp['syngeneic'], pp[c], nbin=20, ax=axs[i], summary_metric='mean')
        axs[i].set_xlabel("# syngeneic - # %s" % c.title())
        i += 1

    axs[-1].legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_histogram_delta_sum_pvalues.png"), dpi=200)

    # plot ECDF (ish)
    p_order = dm_all_p.sum(axis=1).sort_values(ascending=False).index
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for c in comps:
        ax.plot(pp[c].loc[p_order].values.cumsum(), label=c.title())
    ax.set_xlabel('Ranked pathway')
    ax.set_ylabel('Cumulative sum of -log10(p)')
    ax.legend(loc='upper left', frameon=True, facecolor='w', framealpha=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dmr_sum_pvalues_comparisons_ranked_cumul_sum.png"), dpi=200)