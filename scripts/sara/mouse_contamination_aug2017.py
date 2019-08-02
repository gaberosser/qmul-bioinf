"""
WTCHG project P170390: SB had submitted two samples, both commercial cell lines that are supposed to be human
medulloblastoma cultures.

Recall also that we had a mouse contamination with the 1487 xenograft line from Toronto.

One came back normal, the other aligned ~85% with the mouse genome.
This script contains code used to investigate the contamination.
"""
import os
from load_data import rnaseq_data
from plotting import corr
from matplotlib import pyplot as plt
from settings import DATA_DIR
import pandas as pd


# load data

# first, all of the MOUSE samples in the affected run (3 lanes, 20 samples)
# this includes the contaminated sample

indir = os.path.join(DATA_DIR, 'rnaseq', 'wtchg_p170390')

lanedirs = [
    os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
    os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
    os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
]

metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
countdirs = [os.path.join(d, 'mouse', 'star_alignment') for d in lanedirs]
samples1 = (u'CRL3034_METASTATICMB_P4', u'eNSC3med', u'eNSC5med', u'eNSC6med',
           u'eNSC3mouse', u'eNSC5mouse', u'eNSC6mouse', u'mDura3N1mouse',
           u'mDura5N24Amouse', u'mDura6N6mouse', u'mDura3N1human',
           u'mDura5N24Ahuman', u'mDura6N6human')

obj1 = rnaseq_data.all_samples_multilane_loader(
    countdirs, metafiles, source='star', annotate_by='Ensembl Gene ID', samples=samples1
)

# also include the 1487

indir = os.path.join(DATA_DIR, 'rnaseq', 'wtchg_p160704')
lanedirs = [
    os.path.join(indir, '161219_K00198_0151_BHGYHTBBXX'),
    os.path.join(indir, '161222_K00198_0152_AHGYG3BBXX'),
]
metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
countdirs = [os.path.join(d, 'star_alignment_mouse') for d in lanedirs]

obj2 = rnaseq_data.all_samples_multilane_loader(
    countdirs, metafiles, source='star', annotate_by='Ensembl Gene ID',
    samples=('ICb1078', 'ICb1487')
)

obj = rnaseq_data.MultipleBatchLoader([obj1, obj2])

data = obj.data.loc[obj.data.index.str.contains('ENS')]
cpm = data.divide(obj.meta.read_count, axis=1) * 1e6
cpm = cpm.astype(float)
# only keep abundant genes: results in around 12,000 remaining genes
keep = (cpm > 1).sum(axis=1) > 5

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax = corr.plot_correlation_coefficient_array(data.loc[keep], vmin=0.6, ax=ax)
plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
fig.tight_layout()

# 2) What are the correlation coeffs when we look at individual lanes?
# Only include the endogenous NSC in the endogenous medium

indir = os.path.join(DATA_DIR, 'rnaseq', 'wtchg_p170390')
lanedirs = [
    os.path.join(indir, '170727_K00198_0222_AHKWW5BBXX'),
    os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_1'),
    os.path.join(indir, '170731_K00150_0226_AHL2CJBBXX_2')
]
metafiles = [os.path.join(d, 'sources.csv') for d in lanedirs]
countdirs = [os.path.join(d, 'mouse', 'star_alignment') for d in lanedirs]

lane_data = []
lane_cpm = []
lane_samples = []
the_samples = (u'CRL3034_METASTATICMB_P4', u'eNSC3med', u'eNSC5med', u'eNSC6med', u'eNSC3mouse', u'eNSC5mouse')

for i in range(len(countdirs)):
    t = rnaseq_data.StarCountLoader(
        strandedness='r',
        count_dir=countdirs[i],
        meta_fn=metafiles[i],
        samples=the_samples
    )
    this_dat = t.data.loc[t.data.index.str.contains('ENS')]
    this_cpm = this_dat.divide(t.meta.read_count.values, axis=1) * 1e6
    lane_data.append(this_dat)
    lane_cpm.append(this_cpm)
    lane_samples.extend(["%s_%d" % (x, i + 1) for x in the_samples])

# put this in a single DF
per_lane_counts = pd.DataFrame(index=lane_data[0].index, columns=lane_samples)
per_lane_cpm = pd.DataFrame(index=lane_data[0].index, columns=lane_samples)

for i, ct in enumerate(lane_data):
    per_lane_counts.loc[:, ["%s_%d" % (x, i + 1) for x in the_samples]] = ct.values

for i, cp in enumerate(lane_cpm):
    per_lane_cpm.loc[:, ["%s_%d" % (x, i + 1) for x in the_samples]] = cp.values

per_lane_counts = per_lane_counts.loc[:, per_lane_counts.columns.sort_values()]
per_lane_cpm = per_lane_cpm.loc[:, per_lane_cpm.columns.sort_values()]
per_lane_cpm = per_lane_cpm.astype(float)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax = corr.plot_correlation_coefficient_array(per_lane_counts, vmin=0.9, ax=ax)
plt.setp(ax.xaxis.get_ticklabels(), rotation=90)
fig.tight_layout()