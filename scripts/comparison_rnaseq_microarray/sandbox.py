from scripts.comparison_rnaseq_microarray import load_illumina_data, load_rnaseq_data, load_references, comparisons
from microarray.process import aggregate_by_probe_set
import references
import pandas as pd
import numpy as np
from scipy import stats
import collections
from matplotlib import rc, pyplot as plt, gridspec as gridspec
import seaborn as sns
plt.interactive(True)
sns.set_style('white')
# rc('text', usetex=True)


ALL_NORTHCOTT = [
    'WIF1',
    'TNC',
    'GAD1',
    'DKK2',
    'EMX2',
    'ADAMTSL1',
    'NKD1',
    'PDE11A',
    'EPHA7',
    'RUNX2',
    'C20orf103',
    'ABHD12B',
    'PAX3',
    'LGR5',
    'EPHA3',
    'PGM5',
    'TPH1',
    'TNFRSF19',
    'TMEM51',
    'C9orf94',
    'GABRE',
    'TNFRSF11B',
    'FREM2',
    'LEF1',
    'GABRG3',
    'PDLIM3',
    'EYA1',
    'HHIP',
    'ATOH1',
    'SFRP1',
    'PPP2R2C',
    'SOX2',
    'DMRTA1',
    'NDP',
    'GNG3',
    'NDST3',
    'GRIA4',
    'SATB2',
    'CXCR4',
    'CYYR1',
    'SCG5',
    'ALDH1A3',
    'KIF26A',
    'C6orf117',
    'BOC',
    'PRLR',
    'C4orf18',
    'CA14',
    'POU3F2',
    'SEMA6A',
    'IMPG2',
    'GABRA5',
    'EGFL11',
    'NRL',
    'MAB21L2',
    'TTR',
    'NPR3',
    'TBR1',
    'FSTL5',
    'TMEM16B',
    'C5orf23',
    'GNB3',
    'DMD',
    'PCDH21',
    'USH2A',
    'RCVRN',
    'PDE6H',
    'RASGRF2',
    'FOXG1B',
    'SAMD3',
    'TSHZ3',
    'MAK',
    'PPP2R2B',
    'RD3',
    'FAM19A4',
    'KCNA1',
    'EOMES',
    'KHDRBS2',
    'RBM24',
    'UNC5D',
    'OAS1',
    'GRM8',
    'CDH18',
    'LOC138046',
    'SNCAIP',
    'MPP3',
    'CHL1',
    'PEX5L',
    'LMX1A',
    'GPR12',
    'FBXL21',
    'SH3GL3',
    'NID2',
    'LINGO2',
    'PTHLH',
    'CA4',
    'PRL',
    'KCNIP4',
    'NEUROD2',
    'ST18',
    'OTX2',  # not Northcott, requested by SB
]

NORTHCOTT_GENES = (
    ('WNT', ALL_NORTHCOTT[ALL_NORTHCOTT.index('WIF1'):ALL_NORTHCOTT.index('PDLIM3')]),
    ('SHH', ALL_NORTHCOTT[ALL_NORTHCOTT.index('PDLIM3'):ALL_NORTHCOTT.index('IMPG2')]),
    ('Group C', ALL_NORTHCOTT[ALL_NORTHCOTT.index('IMPG2'):ALL_NORTHCOTT.index('KCNA1')]),
    ('Group C', ALL_NORTHCOTT[ALL_NORTHCOTT.index('KCNA1'):]),
)

NANOSTRING_GENES = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1', 'OTX2')),  # OTX2 added by SB
)

# MB_GROUPS = (
#     ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
#     ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
#     ('Group C', ('IMPG2', 'GABRA5', 'EYS', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
#     ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1', 'OTX2')),  # OTX2 added by SB
# )

SAMPLE_GROUPS = (
    ('WNT', ('Pt1140', 'ICb1140-II', 'ICb1140-III', 'Pt1192', 'ICb1192-I', 'ICb1192-III', 'ICb1192-V')),
    ('SSH', ('Pt1338', 'ICb1338-I', 'ICb1338-III', 'ICb984-I', 'ICb984-III', 'ICb984-V')),
    ('Group C', (
        'ICb1197-I',
        'ICb1197-III',
        'Pt1494',
        'ICb1494-I',
        'ICb1494-III',
        'ICb1494-V',
        'Pt1572',
        'ICb1572-I',
        'ICb1572-III',
        'ICb1572-V',
        'Pt1595',
        'ICb1595-I',
        'ICb1595-III',
    )),
    ('Group D', (
        'Pt1078',
        'ICb1078-I',
        'ICb1078-III',
        'ICb1078-V',
        'Pt1299',
        'ICb1299-I',
        'ICb1299-III',
        'ICb1299-IV',
        'Pt1487',
        'ICb1487-I',
        'ICb1487-III',
    )),
)

# lookup Entrez ID of extended Northcott
from references import known_genes

df = known_genes()
all_by_symbol = df.set_index('Approved Symbol').loc[:, ['Entrez Gene ID']].dropna().astype(int)
ncott_id = all_by_symbol.loc[ALL_NORTHCOTT]

matches = df.loc[df.loc[:, 'Approved Symbol'].isin(ALL_NORTHCOTT), 'Approved Symbol'].values
unmatched = ncott_id.index.difference(matches)

# there's a nicer way to do this, but doesn't really matter for this purpose.
ncott_unmatched = {}

for t in unmatched:
    s = df.loc[:, 'Previous Symbols'].dropna().str.contains(t)
    if s.sum() == 1:
        ncott_unmatched[t] = df.loc[s.values, 'Entrez Gene ID'].values[0]
    elif s.sum() == 0:
        ncott_unmatched[t] = np.nan
    else:
        print df.loc[s.values]

unmatches = ncott_id.isnull()

ncott_unmatched = pd.Series(ncott_unmatched)
ncott_id.loc[ncott_unmatched.index, 'Entrez Gene ID'] = ncott_unmatched

unmatched = ncott_id.loc[ncott_id.loc[:, 'Entrez Gene ID'].isnull()].index

ncott_unmatched = {}

for t in unmatched:
    s = df.loc[:, 'Synonyms'].dropna().str.contains(t)
    if s.sum() == 1:
        ncott_unmatched[t] = df.loc[s.values, 'Entrez Gene ID'].values[0]
    elif s.sum() == 0:
        ncott_unmatched[t] = np.nan
    else:
        print df.loc[s.values]

ncott_unmatched = pd.Series(ncott_unmatched)
ncott_id.loc[ncott_unmatched.index, 'Entrez Gene ID'] = ncott_unmatched
