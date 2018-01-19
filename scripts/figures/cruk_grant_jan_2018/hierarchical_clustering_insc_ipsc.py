from load_data import rnaseq_data
from plotting import clustering
from scripts.rnaseq import gtf_reader
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils.output import unique_output_dir
from stats import transformations
import os
from settings import LOCAL_DATA_DIR


def cluster_logdata_with_threshold(data, min_val=None, n=None, mad=None, min_over=2, eps=1e-2, **kwargs):
    func = lambda x: np.log2(data + eps)
    return cluster_data_with_threshold(data, min_val, n, mad=mad, min_over=min_over, transform=func, **kwargs)


def cluster_data_with_threshold(data, min_val=None, n=None, mad=None, min_over=2, transform=None, **kwargs):
    if min_val is not None and min_over is not None:
        idx = (data > min_val).sum(axis=1) > min_over
        data = data.loc[idx]

    if transform is not None:
        data = transform(data)

    if n is not None:
        if mad is None:
            mad = transformations.median_absolute_deviation(data).sort_values(ascending=False)
        else:
            mad = mad.sort_values(ascending=False)
            if len(mad.index.intersection(data.index)) != data.shape[0]:
                raise AttributeError("If a pre-computed MAD is supplied, it must contain all required entries")

        data = data.loc[mad.index[:n]]

    cm = clustering.plot_clustermap(data, cmap='RdBu_r', metric='correlation', **kwargs)
    cm.gs.update(bottom=0.2)
    return cm, mad


if __name__ == "__main__":
    # units = 'estimated_counts'
    units = 'tpm'
    # remove_mt = True
    remove_mt = False
    outdir = unique_output_dir("cruk_figures_hc")
    # load 12 patients iNSC, 4 iPSC
    pids = ['017', '018', '019', '030', '031', '026', '044', '049', '050', '052', '054', '061']
    if units == 'tpm':
        min_val = 1
        min_n = 4
        eps = .01
    elif units == 'estimated_counts':
        min_val = 10
        min_n = 4
        eps = .01

    if remove_mt:
        mt_ensg = set(gtf_reader.get_mitochondrial())

    patient_data = rnaseq_data.load_salmon_by_patient_id(pids, units=units)

    # discard GBM
    patient_data = patient_data.loc[:, ~patient_data.columns.str.contains('GBM')]

    # update index to remove accession version
    idx = patient_data.index.str.replace(r'.[0-9]+$', '')
    patient_data.index = idx

    # now aggregate to gene level and repeat
    # TODO: move to rnaseq module or similar

    fn = os.path.join(LOCAL_DATA_DIR, 'reference_genomes', 'human', 'ensembl', 'GRCh38.p10.release90', 'gene_to_transcript.txt')
    gene_transcript = pd.read_csv(fn, header=0, sep='\t').set_index('Transcript stable ID')

    # shouldn't be necessary, but remove transcripts that have no translation
    to_keep = patient_data.index.intersection(gene_transcript.index)
    if len(to_keep) != patient_data.shape[0]:
        to_drop = patient_data.index.difference(gene_transcript.loc[:, 'Transcript stable ID'])
        print "Discarding %d transcripts that have no associated gene: %s" % (
            len(to_drop), ', '.join(to_drop)
        )
        patient_data = patient_data.loc[to_keep]

    # gene list in same order as data
    genes = gene_transcript.loc[patient_data.index, 'Gene stable ID']

    patient_data_by_gene = patient_data.groupby(genes).sum()

    # discard mitochondrial genes
    if remove_mt:
        idx = ~patient_data_by_gene.index.isin(mt_ensg)
        patient_data_by_gene = patient_data_by_gene.loc[idx]
        # renorm
        if units == 'tpm':
            patient_data_by_gene = patient_data_by_gene.divide(patient_data_by_gene.sum(), axis=1) * 1e6

    # bring in reference data
    # IDs (if req), lab (appears in label), loader
    ref_dats = [
        (None, 'Caren et al.', rnaseq_data.pollard_salmon(units=units),),
        (None, 'Yang et al.', rnaseq_data.gse80732_salmon(units=units),),
        (None, 'Shahbazi et al.', rnaseq_data.gse64882_salmon(units=units),),
        (None, 'Li et al.', rnaseq_data.gse84166_salmon(units=units),),
        (None, 'Duan et al.', rnaseq_data.gse61794_salmon(units=units),),
    ]

    ref_arr = []
    ref_labels = []
    ref_cols = []
    batches = []

    for i, r in enumerate(ref_dats):
        the_dat = r[2]
        if r[0] is not None:
            the_cols = r[0]
        else:
            the_cols = the_dat.columns
        ref_cols.extend(the_cols)
        ref_labels.extend(the_dat.columns)
        ref_arr.append(the_dat)
        batches.extend([r[1]] * the_dat.shape[1])

    ref = pd.concat(ref_arr, axis=1)
    ref.columns = ref_cols
    batches = pd.Series(batches, index=ref_cols)
    labels = pd.Series(ref_labels, ref_cols)
    ref.index = ref.index.str.replace(r'.[0-9]+$', '')

    # discard immortalised cell line
    # discard fibroblasts (careful, need to be selective here)

    to_discard = [
        'Fetal ctx', 'tumor', 'ctx endo', 'whole cortex', 'myeloid', 'LM-NSC', r'^fibroblast', 'NHF1-hTERT fibroblasts',
        'foreskin fibroblast', 'lung fibroblast'
    ]
    for td in to_discard:
        the_idx = ~ref.columns.str.contains(td)
        ref = ref.loc[:, the_idx]
        batches = batches.loc[the_idx]
        labels = labels.loc[the_idx]

    # Take mean over various replicates (for simplicity)
    to_aggr = [
        ('Pollard NSC', 'Fetal NSC'),
        (r'H1-[12]', 'H1 PSC'),
        (r'H1-EPS', 'H1 EPSC'),
        (r'ES1-EPS', 'ES1 EPSC'),
        (r'H9_NSC', 'H9 NSC'),
        (r'INSC fibroblast', 'Fibroblast-derived iNSC'),
        (r'NSC008', 'Fetal NSC'),
        (r'fetal NSC rep', 'Fetal NSC'),
    ]

    for srch, repl in to_aggr:
        idx = ref.columns.str.contains(srch)
        if idx.sum() == 0:
            continue

        new_col = ref.loc[:, idx].mean(axis=1)
        new_batch = batches.loc[idx].iloc[0]
        new_label = repl

        ref = ref.loc[:, ~idx]
        ref.insert(ref.shape[1], repl, new_col, allow_duplicates=True)
        batches = batches.loc[~idx]
        batches = batches.append(pd.Series([new_batch], index=[repl]))
        labels = labels.loc[~idx]
        labels = labels.append(pd.Series([new_label], index=[repl]))

    # append study name to labels
    new_labels = ["%s (%s)" % (t) for t in zip(labels.index, batches.values)]

    ref.columns = new_labels
    batches.index = new_labels

    ref_by_gene = ref.groupby(genes).sum()

    # now let's try clustering everything together
    abg = pd.concat((patient_data_by_gene, ref_by_gene), axis=1)

    # discard mitochondrial genes
    if remove_mt:
        idx = ~abg.index.isin(mt_ensg)
        abg = abg.loc[idx]
        # renorm
        if units == 'tpm':
            abg = abg.divide(abg.sum(), axis=1) * 1e6

    # discard genes expressed at low values
    idx = (abg > min_val).sum(axis=1) > min_n
    abg = abg.loc[idx]

    abg_log = np.log2(abg + eps)
    amad_log = transformations.median_absolute_deviation(abg_log).sort_values(ascending=False)

    if units == 'estimated_counts':
        # optionally could normalise here?
        pass

    row_colours_all = pd.DataFrame('gray', index=abg.columns, columns=[''])
    row_colours_all.loc[row_colours_all.index.str.contains(r'NSC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'NPC')] = 'blue'
    row_colours_all.loc[row_colours_all.index.str.contains(r'GIBCO')] = '#96daff'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Fibroblast')] = '#fff89e'
    row_colours_all.loc[row_colours_all.index.str.contains(r'Fetal')] = 'yellow'
    row_colours_all.loc[row_colours_all.index.str.contains('ES1')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains('PSC')] = '#ff7777'
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_NSC')] = '#7fc97f'  # green
    row_colours_all.loc[row_colours_all.index.str.contains(r'DURA[0-9]*_IPSC')] = '#fdc086'  # orange

    n_gene = 3000

    fname = "hier_clust_by_gene_log_corr_top%d_by_mad.{ext}" % n_gene

    cm, mad_all = cluster_logdata_with_threshold(abg, n=n_gene, eps=eps, col_colors=row_colours_all)
    cm.gs.update(bottom=0.3)
    cm.savefig(os.path.join(outdir, fname.format(ext='png')), dpi=300)
    cm.savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)

    fname = "hier_clust_dendrogram_log_corr_top%d_by_mad.{ext}" % n_gene
    d = clustering.dendrogram_with_colours(
        abg_log.loc[amad_log.index[:n_gene]],
        row_colours_all,
        fig_kws={'figsize': (5.5, 10)},
        vertical=False
    )
    d['fig'].savefig(os.path.join(outdir, fname.format(ext='png')), dpi=300)
    d['fig'].savefig(os.path.join(outdir, fname.format(ext='tiff')), dpi=200)
