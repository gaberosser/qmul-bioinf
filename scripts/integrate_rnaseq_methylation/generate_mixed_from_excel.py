import pandas as pd


infile = '/home/gabriel/Dropbox/research/qmul/results/gbm_insc_paired_sample_de/all/n_equals_2/gene_lists/individual_gene_lists_de_dmr.xlsx'
dat = pd.read_excel(infile, None)

mixed = {}
cols2 = ['gene',
 'de_logfc',
 'de_logcpm',
 'de_padj',
 'dmr_cid',
 'dmr_chr',
 'dmr_median1',
 'dmr_median2',
 'dmr_median_delta',
 'dmr_padj',
 'dmr_class_gene',
 'dmr_class_island',
 'dmr_class_tss',
 'de_direction',
 'dmr_direction',
 'de_dmr_concordant']

cols = pd.MultiIndex.from_tuples(
    [('pair', t) for t in cols2] + [('ref', t) for t in cols2],
    names=['comparison', 'fields']
)

pids = set()
for k in dat:
    pids.add(k[3:6])
pids = list(pids)

mixed = {}
for pid in pids:
    the_dat = dat["GBM%s_pair_and_ref" % pid].iloc[2:]
    dat_pair = the_dat.iloc[:, 1:17]
    dat_pair.columns = cols2
    dat_ref = the_dat.iloc[:, 17:]
    dat_ref.columns = cols2

    dat_pair_id = dat_pair.loc[:, ['gene', 'dmr_cid']].apply(tuple, axis=1)
    dat_ref_id = dat_ref.loc[:, ['gene', 'dmr_cid']].apply(tuple, axis=1)

    dat_pair.index = dat_pair_id
    dat_ref.index = dat_ref_id

    dat_ref = dat_ref.loc[dat_pair.index]

    this_idx = dat_pair.loc[:, 'de_dmr_concordant'] != dat_ref.loc[:, 'de_dmr_concordant']

    this_mixed = pd.concat((dat_pair.loc[this_idx], dat_ref.loc[this_idx]), axis=1)
    this_mixed.columns = cols
    this_mixed.index = range(this_mixed.shape[0])

    mixed[pid] = this_mixed

out_fn = 'de_dmr_mixed_concordance.xlsx'
xl_writer = pd.ExcelWriter(out_fn)

# sort the keys for a more sensible order
keys = sorted(mixed.keys())
for k in keys:
    bl = mixed[k]
    bl.to_excel(xl_writer, k, index=True)
xl_writer.save()