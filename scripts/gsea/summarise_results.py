import pandas as pd
import os
from glob import glob


if __name__ == '__main__':
    base_dir = '/home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gsea/tpm_salmon/geneset_permutation/'

    pids = ['017_v2']  # TODO: include all pids
    fdr = 0.05

    for pid in pids:
        pathways = {}
        indir = os.path.join(base_dir, pid)
        the_report = glob(os.path.join(indir, "gsea_report_for_GBM*.xls"))[0]
        report = pd.read_csv(the_report, sep='\t', header=0, index_col=None, usecols=[0, 3, 4, 5, 6, 7, 8])
        report.columns = [
            'pathway',
            'n_gene',
            'es',
            'nes',
            'nom_pval',
            'fdr',
            'fwer',
        ]
        # filter by FDR
        report = report.loc[report.fdr <= fdr]

        for p in report.pathway.values:
            infile = os.path.join(indir, "%s.xls" % p)
            if not os.path.isfile(infile):
                print "Unable to find expected file %s. Skipping." % infile
            else:
                the_path = pd.read_csv(infile, sep='\t', header=0, index_col=None, usecols=[1, 8])
                pathways[p] = the_path

        all_genes = reduce(lambda x, y: set(x).union(y), (v.PROBE for v in pathways.values()))
        all_genes = sorted(all_genes)

        # create the combined result
        res = pd.DataFrame(0, index=pathways.keys(), columns=all_genes)
        for p, v in pathways.items():
            enriched = v.loc[:, 'CORE ENRICHMENT'] == 'Yes'
            not_enriched = ~enriched
            res.loc[p, v.PROBE[enriched].values] = 2
            res.loc[p, v.PROBE[not_enriched].values] = 1

        report.set_index('pathway', inplace=True)
        res.insert(0, 'fdr', report.loc[res.index, 'fdr'])
        ss = res.fdr.sort_values()
        res = res.loc[ss.index]
        res.to_csv(os.path.join(indir, "%s_summary.csv" % pid))