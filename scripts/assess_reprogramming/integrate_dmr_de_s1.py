from scripts.assess_reprogramming import analyse_residual_de_novo_results as ardnr
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts
from utils import output, log
from rnaseq import loader
from settings import OUTPUT_DIR

import os
import pickle
import re
import pandas as pd
import statsmodels.formula.api as sm

logger = log.get_console_logger()


if __name__ == '__main__':
    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'de')

    rnaseq_obj = loader.load_by_patient(consts.PIDS)
    rnaseq_obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES)

    dat_s1 = rnaseq_obj.data
    meta_s1 = rnaseq_obj.meta.loc[dat_s1.columns]

    the_hash = tsgd.de_results_hash(meta_s1.index.tolist(), consts.DE_PARAMS)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_full_s1 = pickle.load(f)
    else:
        raise NotImplementedError("We require a precomputed DE results file (expected this at %s)" % fn)

    basedir = os.path.join(OUTPUT_DIR, "assess_reprog_alt1_apocrita")
    indir = os.path.join(basedir, "results")
    names = [t.replace('.csv', '') for t in os.listdir(indir) if '.csv' in t]

    fdr = 0.01

    colour_by_direction = {
        'Hypermethylated': '#e09191',  # red
        'Hypomethylated': '#91e097',  # green
    }

    colour_by_pid = {
        '019': '#7fc97f',
        '030': '#beaed4',
        '031': '#fdc086',
        '050': '#ffff99',
        '054': '#386cb0',
    }

    # split into comparison types
    ipsc_vs_fb_names = [t for t in names if re.match(r'iPSC[^-]*-FB.*', t)]
    ipsc_vs_esc_names = [t for t in names if re.match(r'iPSC[^-]*-ESC.*', t)]
    insc_vs_ipsc_names = [t for t in names if re.match(r'iNSC[^-]*-iPSC.*', t)]
    esc_vs_esc_names = [t for t in names if re.match(r'ESC[^-]*-ESC.*', t)]

    all_types = [
        ipsc_vs_esc_names,
        ipsc_vs_fb_names,
        insc_vs_ipsc_names,
        esc_vs_esc_names
    ]

    if sum([len(t) for t in all_types]) != len(names):
        logger.warn("Some comparisons not covered")

    # now load the actual results
    ipsc_ours_regex = r'iPSC0[135].*'
    fb_ours_regex = r'FB0[135]*'
    ipsc_e6194_regex = r'iPSCHEL[^-]*'
    fb_e6194_regex = r'FB27_HFF'
    ipsc_hipsci_regex = r'iPSCHPS[^-]*'
    ipsc_banov_regex = r'iPSCNA[0-9].*'
    esc_encode_regex = r'ESCH7 hESC'
    esc_e6194_regex = r'ESCH9_p50'
    insc_ours_regex = r'iNSC0[135].*'

    ipsc_regex = {
        'ours': ipsc_ours_regex,
        'e6194': ipsc_e6194_regex,
        'hipsci': ipsc_hipsci_regex,
        'banov': ipsc_banov_regex
    }

    fb_regex = {
        'ours': fb_ours_regex,
        'e6194': fb_e6194_regex
    }

    esc_regex = {
        'encode': esc_encode_regex,
        'e6194': esc_e6194_regex
    }

    insc_regex = {
        'ours': insc_ours_regex
    }

    ipsc_vs_fb = {}
    for k1, v1 in ipsc_regex.items():
        for k2, v2 in fb_regex.items():
            the_regex = re.compile("%s-%s" % (v1, v2))
            the_names = [t for t in ipsc_vs_fb_names if re.match(the_regex, t)]

            ipsc_vs_fb.setdefault(k1, {})[k2] = dict([
                (t, pd.read_csv(os.path.join(indir, "%s.csv" % t), header=0, index_col=0)) for t in the_names
            ])

    ipsc_vs_esc = {}
    for k1, v1 in ipsc_regex.items():
        for k2, v2 in esc_regex.items():
            the_regex = re.compile("%s-%s" % (v1, v2))
            the_names = [t for t in ipsc_vs_esc_names if re.match(the_regex, t)]

            ipsc_vs_esc.setdefault(k1, {})[k2] = dict([
                (t, pd.read_csv(os.path.join(indir, "%s.csv" % t), header=0, index_col=0)) for t in the_names
            ])

    insc_vs_ipsc = {}
    for k1, v1 in insc_regex.items():
        for k2, v2 in ipsc_regex.items():
            the_regex = re.compile("%s-%s" % (v1, v2))
            the_names = [t for t in insc_vs_ipsc_names if re.match(the_regex, t)]

            insc_vs_ipsc.setdefault(k1, {})[k2] = dict([
                (t, pd.read_csv(os.path.join(indir, "%s.csv" % t), header=0, index_col=0)) for t in the_names
            ])

    esc_vs_esc = {}
    for k1, v1 in esc_regex.items():
        for k2, v2 in esc_regex.items():
            the_regex = re.compile("%s-%s" % (v1, v2))
            the_names = [t for t in esc_vs_esc_names if re.match(the_regex, t)]

            esc_vs_esc.setdefault(k1, {})[k2] = dict([
                (t, pd.read_csv(os.path.join(indir, "%s.csv" % t), header=0, index_col=0)) for t in the_names
            ])

    # get lists of sample names
    ipsc_samples = {}
    fb_samples = {}
    esc_samples = {}
    insc_samples = {}

    for i in ipsc_vs_fb:
        for j in ipsc_vs_fb[i]:
            all_ipsc = sorted(set([t.split('-')[0] for t in ipsc_vs_fb[i][j]]))
            all_fb = sorted(set([t.split('-')[1] for t in ipsc_vs_fb[i][j]]))
            ipsc_samples[i] = all_ipsc
            fb_samples[j] = all_fb

    for i in ipsc_vs_esc:
        for j in ipsc_vs_esc[i]:
            all_esc = sorted(set([t.split('-')[1] for t in ipsc_vs_esc[i][j]]))
            esc_samples[j] = all_esc

    for i in insc_vs_ipsc:
        all_insc = sorted(set([t.split('-')[0] for t in insc_vs_ipsc[i]['ours']]))
        insc_samples[i] = all_insc

    # Our DMRs: only compare with parental lines
    ipsc_esc_core = ardnr.core_dmrs(ipsc_vs_esc, fdr=fdr)

    ipsc_fb_ours = {}
    for k, df in ipsc_vs_fb['ours']['ours'].items():
        ii, jj = k.split('-')
        if ii.replace('iPSC', '') == jj.replace('FB', ''):
            ipsc_fb_ours[ii] = df.loc[df.padj < fdr]

    # UPDATE: no longer excluding DMRs that are also DMRs in the ESC-ESC comparison
    # In practice, this makes little difference
    # core_dmr_classified_ours = classify_dmrs_residual_denovo(ipsc_esc_core['ours'], ipsc_fb_ours, exclude=esc_esc_dmrs)
    core_dmr_classified_ours = ardnr.classify_dmrs_residual_denovo(ipsc_esc_core['ours'], ipsc_fb_ours, exclude=None)
    core_dmr_classified_count_ours = dict(
        [(k, v.classification.value_counts()) for k, v in core_dmr_classified_ours.items()]
    )

    # export to CSV format (one file per patient)
    # also generate a list of genes linked to DMRs in residual hypomethylated regions
    residual_hypo_genes = {}
    for k, v in core_dmr_classified_ours.items():
        v.to_csv(os.path.join(outdir, "%s_classified_dmrs.csv" % k))
        genes = sorted(set(
            v.loc[v.classification == 'hypo_residual', 'genes'].apply(ardnr.make_tuple).sum()
        ))
        residual_hypo_genes[k] = genes

    # Other matching comparisons
    ipsc_fb_e6194 = {}
    for k, df in ipsc_vs_fb['e6194']['e6194'].items():
        ii, jj = k.split('-')
        if re.search(r'HEL(140|141)', ii):
            ipsc_fb_e6194[ii] = df.loc[df.padj < fdr]

    # core_dmr_classified_e6194 = classify_dmrs_residual_denovo(ipsc_esc_core['e6194'], ipsc_fb_e6194, exclude=esc_esc_dmrs)
    core_dmr_classified_e6194 = ardnr.classify_dmrs_residual_denovo(ipsc_esc_core['e6194'], ipsc_fb_e6194, exclude=None)
    core_dmr_classified_count_e6194 = dict(
        [(k, v.classification.value_counts()) for k, v in core_dmr_classified_e6194.items()]
    )


    ## TODO: fix this mess and integrate DE/DMR.

    # de_logfc = {}
    # dmr_median_delta = {}
    #
    # for pid in consts.PIDS:
    #     fn = os.path.join(dmr_indir, "iPSC%s_classified_dmrs.csv" % pid)
    #     if (pid in dmrs_classified) and ("iPSC_%s_ours" % pid in ipsc_esc_fb):
    #         this_dmr = dmrs_classified[pid]
    #         dmr_genes = this_dmr.loc[:, 'genes'].values
    #         dmr_genes = setops.reduce_union(*dmr_genes) if len(dmr_genes) else []
    #         this_de_sign = ipsc_esc_fb["iPSC_%s_ours" % pid]
    #         this_de_full = {}
    #         for r in ref_labels:
    #             tt = de_res_full[("iPSC_%s_ours" % pid, "ESC_%s" % r)]
    #             tt = tt.loc[tt['Gene Symbol'].isin(dmr_genes)]
    #             this_de_full[r] = tt
    #
    #         common_ix = sorted(setops.reduce_intersection(*[t.index for t in this_de_full.values()]))
    #         common_gs = this_de_full.values()[0].loc[common_ix, 'Gene Symbol']
    #
    #         this_dmr_dm = {}
    #
    #         for g in common_gs:
    #             this_ix = this_dmr.index[this_dmr.genes.apply(lambda x: g in x)]
    #             this_dmr_dm[g] = this_dmr.loc[this_ix, this_dmr.columns.str.contains('median_delta')].mean(axis=0)
    #         dmr_median_delta[pid] = pd.DataFrame(this_dmr_dm).transpose().sort_index()
    #
    #         this_logfc = pd.concat(
    #             (this_de_full[r].loc[common_ix, 'logFC'] for r in ref_labels),
    #             axis=1
    #         )
    #         this_logfc.columns = ref_labels
    #         this_logfc.index = common_gs
    #
    #         de_logfc[pid] = this_logfc.sort_index()
    #
    #
    # for pid in consts.PIDS:
    #     if pid in de_logfc:
    #         fig = plt.figure(figsize=(5, 4))
    #         ax = fig.add_subplot(111)
    #         this_de = de_logfc[pid].mean(axis=1)
    #         this_dm = dmr_median_delta[pid].mean(axis=1)
    #         ax.scatter(this_dm, this_de, edgecolor='k', facecolor='royalblue', zorder=2)
    #         ax.errorbar(
    #             this_dm,
    #             this_de,
    #             xerr=dmr_median_delta[pid].diff(axis=1).dropna(axis=1).values / 2.,
    #             yerr=de_logfc[pid].diff(axis=1).dropna(axis=1).values / 2.,
    #             ls='none',
    #             color='royalblue',
    #             zorder=1
    #         )
    #         ws = pd.DataFrame({'x': this_dm, 'y': this_de})
    #         xerr = dmr_median_delta[pid].diff(axis=1).dropna(axis=1).values
    #         yerr = de_logfc[pid].diff(axis=1).dropna(axis=1).values
    #         weights = pd.Series(((xerr ** 2 + yerr ** 2) ** -.5).flatten(), index=ws.index)
    #         wls_fit = sm.wls('y ~ x', data=ws, weights=weights).fit()
    #         print "Patient %s" % pid
    #         print wls_fit.summary()
    #         ax.plot(ws.x, wls_fit.predict(), color='darkgray', linewidth=1.5)
    #         ax.axhline(0, linestyle='--', color='k', linewidth=1.)
    #         ax.axvline(0, linestyle='--', color='k', linewidth=1.)
    #         ax.set_xlabel(r'Methylation median $\Delta M$')
    #         ax.set_ylabel(r'Gene expression $\log_2$ fold change')
    #         fig.tight_layout()
    #         fig.savefig(os.path.join(outdir, "de_logfc_vs_dmr_delta_m_%s.png" % pid), dpi=200)
    #         fig.savefig(os.path.join(outdir, "de_logfc_vs_dmr_delta_m_%s.tiff" % pid), dpi=200)
    #
    # # test: are these DMRs significantly more concordant with DE than randomly chosen ones??
    # # statistic to use: rsq of weighted linear regression
    # # This will replace the statistics reported by WLS (to some extent)
    # ## TODO!
    #
    # n_perm = 1000
    # rsq_true = {}
    # rsq_permute = {}
    #
    # for pid in de_logfc:
    #     this_de = de_logfc[pid].mean(axis=1)
    #     this_dm = dmr_median_delta[pid].mean(axis=1)
    #     ws = pd.DataFrame({'x': this_dm, 'y': this_de})
    #     xerr = dmr_median_delta[pid].diff(axis=1).dropna(axis=1).values
    #     yerr = de_logfc[pid].diff(axis=1).dropna(axis=1).values
    #     weights = pd.Series(((xerr ** 2 + yerr ** 2) ** -.5).flatten(), index=ws.index)
    #     wls_fit = sm.wls('y ~ x', data=ws, weights=weights).fit()
    #     rsq_true[pid] = rsq_true
    #
    #     # now repeat for randomly-selected genes