import os
from settings import HGIC_LOCAL_DIR
import pandas as pd
from utils import output, setops, excel
from plotting import common
import consts
from matplotlib import pyplot as plt, patches, collections, gridspec
import seaborn as sns
import numpy as np


def pathway_involvement_heatmap_by_p(
        p_dat,
        n_set,
        pathway_order,
        pids,
        comparison_dict,
        orientation='vertical',
        vmin=None,
        vmax=None,
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
        cmap='Blues',
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


if __name__ == '__main__':
    # set a minimum pval for pathways to be used
    alpha = 0.005
    plogalpha = -np.log10(alpha)
    # more lenient pval threshold for considering pathways as relevant
    alpha_relevant = 0.05
    plogalpha_relevant = -np.log10(alpha_relevant)

    indir = os.path.join(HGIC_LOCAL_DIR, 'current/core_pipeline/rnaseq/merged_s1_s2/ipa/pathways')
    outdir = output.unique_output_dir()

    # keys are the term used in the filename, values are those used in the columns
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

    # first, load from raw data and combine into a single export file
    # format: Excel, wideform
    file_patt = 'de_s2_{pid}_{cmp}.txt'
    to_export = {}
    col_order = []
    for pid in pids:
        for c in comps:
            fn = os.path.join(indir, file_patt.format(pid=pid, cmp=c))
            this = pd.read_csv(fn, sep='\t', skiprows=2, header=0, index_col=0)
            this.columns = ['-logp', 'ratio', 'z', 'genes']
            # replace genes column with n genes (to avoid overcomplicating it)
            this.insert(3, 'n_gene', this.genes.str.split(',').apply(len))
            # resolve encoding in index
            this.index = [x.decode('utf-8') for x in this.index]
            # restrict to relevant pathways
            rele_ix = this.index[this['-logp'] >= plogalpha_relevant]
            this = this.loc[rele_ix]
            to_export["%s_%s" % (pid, c)] = this
            col_order.append("%s_%s" % (pid, c))

    # wideform version of this (i.e. 30 blocks)
    # we can't use the Venn approach here, but we don't need to
    all_pathways = sorted(setops.reduce_union(*[t.index for t in to_export.values()]))
    export_wideform = pd.DataFrame(index=all_pathways)
    member_cols = []
    for pid in pids:
        for c in comps:
            k = "%s_%s" % (pid, c)
            this = to_export[k]
            sign_ix = this.index[this['-logp'] >= plogalpha]
            this_yn = pd.Series('N', index=all_pathways)
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
                    this.reindex(all_pathways)[col]
                )

    # add n gene in pathway as single const column
    rr = export_wideform.loc[:, export_wideform.columns.str.contains('ratio')]
    ng = export_wideform.loc[:, export_wideform.columns.str.contains('n_gene')]
    n_gene_tot = (ng.astype(float).values / rr.astype(float)).mean(axis=1).round().astype(int)
    export_wideform.insert(0, 'n_gene_in_pathway', n_gene_tot)

    export_wideform.to_excel(os.path.join(outdir, "full_de_ipa_results.xlsx"))

    # first pass: load data and obtain list of pathways considered 'significant' (based on alpha)
    pathways_to_keep = set()

    ipa_res = {}
    ipa_res_full = {}
    for c, kc in comps.items():
        p = pd.read_csv(
            os.path.join(indir, 'full_de_%s_logp.txt' % c),
            sep='\t',
            skiprows=2,
            index_col=0,
            header=0,
            usecols=range(11)
        )
        z = pd.read_csv(
            os.path.join(indir, 'full_de_%s_z.txt' % c),
            sep='\t',
            skiprows=2,
            index_col=0,
            header=0,
            usecols=range(11)
        ).reindex(p.index)
        p.columns = p.columns.str.replace('_logFC', '')
        z.columns = z.columns.str.replace('_logFC', '')
        for pid in pids:
            t = pd.DataFrame(index=p.index, columns=['-logp', 'z'])
            t.loc[:, '-logp'] = p["%s_%s" % (pid, kc)]
            t.loc[:, 'z'] = z["%s_%s" % (pid, kc)]
            # resolve encoding in index
            t.index = [x.decode('utf-8') for x in t.index]

            # store relevant results
            ipa_res_full['%s_%s' % (pid, c)] = t.loc[t['-logp'] > -np.log10(alpha_relevant)]
            # store significant results
            ipa_res['%s_%s' % (pid, c)] = t.loc[t['-logp'] > -np.log10(alpha)]

            pathways_to_keep.update(ipa_res['%s_%s' % (pid, c)].index)

    excel.pandas_to_excel(ipa_res, os.path.join(outdir, "full_de_ipa_results_significant_separated.xlsx"))

    pathways_to_keep = sorted(pathways_to_keep)

    # use this list to export a second wideform Excel file with the top list of pathways
    p_for_export = pd.DataFrame(index=pathways_to_keep, columns=sorted(ipa_res.keys()))
    z_for_export = pd.DataFrame(index=pathways_to_keep, columns=sorted(ipa_res.keys()))
    s_for_export = pd.DataFrame('N', index=pathways_to_keep, columns=sorted(ipa_res.keys()))
    for k, v in ipa_res_full.items():
        p_for_export.loc[pathways_to_keep, k] = v.reindex(pathways_to_keep).loc[:, '-logp']
        z_for_export.loc[pathways_to_keep, k] = v.reindex(pathways_to_keep).loc[:, 'z']
        s_for_export.loc[v.reindex(pathways_to_keep).loc[:, '-logp'] > -np.log10(alpha), k] = 'Y'
    p_for_export.columns = ["%s_-logp" % t for t in p_for_export.columns]
    z_for_export.columns = ["%s_z" % t for t in z_for_export.columns]
    for_export = pd.concat((p_for_export, z_for_export, s_for_export), axis=1).sort_index(axis=1)

    for_export.to_excel(os.path.join(outdir, "full_de_ipa_results_significant.xlsx"))


    # second pass: obtain full data from each comparison, providing the pathway is 'relevant' (based on alpha_relevant)
    all_p = {}
    all_z = {}
    all_in = {}

    for i, pid in enumerate(pids):
        this_data = []
        this_p = []
        this_z = []
        for c in sorted(comps.keys()):
            the_dat = ipa_res['%s_%s' % (pid, c)].reindex(pathways_to_keep)
            the_dat = the_dat.loc[the_dat['-logp'] > -np.log10(alpha_relevant)]

            all_in["%s_%s" % (pid, c)] = the_dat['-logp'] > np.log10(alpha)
            all_p["%s_%s" % (pid, c)] = the_dat['-logp']
            all_z["%s_%s" % (pid, c)] = the_dat['z']

    all_in = pd.DataFrame(all_in).fillna(False)
    all_z = pd.DataFrame(all_z)
    all_p = pd.DataFrame(all_p)

    # number syngen. only, ref. only and intersection
    n_set = pd.DataFrame(0, index=all_p.index, columns=['Syngen. only', 'Ref. only', 'Intersect.'], dtype=int)
    so = dict([(pw, []) for pw in pathways_to_keep])
    ro = dict([(pw, []) for pw in pathways_to_keep])
    inters = dict([(pw, []) for pw in pathways_to_keep])
    for pid in pids:
        s = all_in.index[all_in.loc[:, "%s_syngeneic" % pid]]
        r = all_in.index[all_in.loc[:, ["%s_%s" % (pid, t) for t in comps.keys()[1:]]].any(axis=1)]
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
        index=pathways_to_keep,
        columns=['n_syngen_only', 'syngen_only_pids', 'n_ref_only', 'ref_only_pids', 'n_intersect', 'intersect_pids'],
        dtype=object
    )
    for pw in pathways_to_keep:
        at_a_glance.loc[pw, 'n_syngen_only'] = len(so[pw])
        at_a_glance.loc[pw, 'syngen_only_pids'] = ';'.join(so[pw])
        at_a_glance.loc[pw, 'n_ref_only'] = len(ro[pw])
        at_a_glance.loc[pw, 'ref_only_pids'] = ';'.join(ro[pw])
        at_a_glance.loc[pw, 'n_intersect'] = len(inters[pw])
        at_a_glance.loc[pw, 'intersect_pids'] = ';'.join(inters[pw])

    at_a_glance.to_excel(os.path.join(outdir, "ipa_results_patients_by_s2_category.xlsx"))

    z_min = all_z.min().min()
    z_max = all_z.max().max()

    # plot 1) P values, ordered by sum of -log(P)
    p_order = all_p.sum(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        all_p,
        n_set,
        p_order,
        pids,
        comparison_names
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_sum_logp_de.tiff"), dpi=200)

    # plot 2) P values, ordered by mean of -log(P)
    p_order = all_p.mean(axis=1).sort_values(ascending=False).index
    plot_dict = pathway_involvement_heatmap_by_p(
        all_p,
        n_set,
        p_order,
        pids,
        comparison_names
    )
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.png"), dpi=200)
    plot_dict['figure'].savefig(os.path.join(outdir, "heatmap_all_pathways_order_mean_logp_de.tiff"), dpi=200)

    # TODO: use this to make a Z value function, then delete?
    if False:

        # again but vertical and with labels
        plot = 'p'


        fig = plt.figure(figsize=(9., 11.5))
        gs_kw = dict(
            left=0.43,
            right=0.93,
            top=0.998,
            bottom=0.1,
            wspace=0.1,
            hspace=0.1,
        )

        # set up axis grid
        gs = gridspec.GridSpec(
            ncols=len(pids) + 2,
            nrows=5,
            width_ratios=[1] * len(pids) + [1, .5],
            **gs_kw
        )
        cax = fig.add_subplot(gs[1:-1, -1])
        axs = []

        for i, pid in enumerate(pids):
            ax = fig.add_subplot(gs[:, i])
            axs.append(ax)
            this_p = []
            this_z = []
            for j, c in enumerate(sorted(comps.keys())):
                # embed this set of data in full sorted list
                the_z = all_z.loc[:, '%s_%s' % (pid, c)].reindex(p_order)
                the_p = all_p.loc[:, '%s_%s' % (pid, c)].reindex(p_order)
                this_z.append(the_z)
                this_p.append(the_p)

            this_z = pd.concat(this_z, axis=1).iloc[:, ::-1]
            this_p = pd.concat(this_p, axis=1).iloc[:, ::-1]

            if plot == 'z':

                # Z not null, P not null: pathway is enriched and has direction
                sns.heatmap(
                    this_z,
                    mask=this_p.isnull() | this_z.isnull(),
                    cmap='RdBu_r',
                    linewidths=.2,
                    linecolor='w',
                    vmin=z_min,
                    vmax=z_max,
                    yticklabels=False,
                    cbar=i == 0,
                    cbar_ax=cax,
                    cbar_kws={"orientation": "vertical"},
                    ax=ax,
                )
                # Z null, P not null: no direction information, but pathway is enriched
                tmp = this_z.fillna(1.)
                tmp[tmp != 1] = 0.
                sns.heatmap(
                    tmp,
                    mask=(~this_z.isnull()) | this_p.isnull(),
                    cmap='binary',
                    linewidths=.2,
                    linecolor='w',
                    yticklabels=False,
                    cbar=False,
                    ax=ax
                )
                cax.set_ylabel('Z score')
            else:
                sns.heatmap(
                    this_p,
                    mask=this_p.isnull(),
                    cmap='YlOrRd',
                    linewidths=.2,
                    linecolor='w',
                    # vmin=z_min,
                    # vmax=z_max,
                    yticklabels=False,
                    cbar=i == 0,
                    cbar_ax=cax,
                    cbar_kws={"orientation": "vertical"},
                    ax=ax,
                )
                cax.set_ylabel('$-\log(p)$')
            plt.setp(ax.xaxis.get_ticklabels(), rotation=0)
            ax.xaxis.set_ticklabels(['Syngen.', 'H9', 'Gibco'], rotation=90.)
            ax.set_xlabel(pid)
            ax.set_facecolor('0.6')

        eax = fig.add_subplot(gs[:, len(pids)])
        sns.heatmap(
            n_set,
            mask=n_set==0,
            cmap='Blues',
            # yticklabels=False,
            cbar=False,
            ax=eax,
            annot=True,
            fmt="d",
            annot_kws={"size": 6}
        )
        eax.yaxis.set_ticks([])
        plt.setp(eax.xaxis.get_ticklabels(), rotation=90)

        # colorbar outline
        cbar = axs[0].collections[0].colorbar
        cbar.outline.set_linewidth(1.)
        cbar.outline.set_edgecolor('k')

        axs[0].set_yticks(0.5 + np.arange(len(p_order)))
        axs[0].set_yticklabels(p_order[::-1], rotation=0, fontsize=7)

        fig.savefig(os.path.join(outdir, "heatmap_all_pathways_de.png"), dpi=200)
        fig.savefig(os.path.join(outdir, "heatmap_all_pathways_de.tiff"), dpi=200)