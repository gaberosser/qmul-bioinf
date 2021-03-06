import collections


NORTHCOTT_GENEID_MAP = dict([
    ('WIF1', 11197),
    ('TNC', 3371),
    ('GAD1', 2571),
    ('DKK2', 27123),
    ('EMX2', 2018),
    ('ADAMTSL1', 92949),
    ('NKD1', 85407),
    ('PDE11A', 50940),
    ('EPHA7', 2045),
    ('RUNX2', 860),
    ('C20orf103', 24141),
    ('ABHD12B', 145447),
    ('PAX3', 5077),
    ('LGR5', 8549),
    ('EPHA3', 2042),
    ('PGM5', 5239),
    ('TPH1', 7166),
    ('TNFRSF19', 55504),
    ('TMEM51', 55092),
    # ('C9orf94', 92949),  # this is a synonym of ADAMTSL1
    ('GABRE', 2564),
    ('TNFRSF11B', 4982),
    ('FREM2', 341640),
    ('LEF1', 51176),
    ('GABRG3', 2567),
    ('PDLIM3', 27295),
    ('EYA1', 2138),
    ('HHIP', 64399),
    ('ATOH1', 474),
    ('SFRP1', 6422),
    ('PPP2R2C', 5522),
    ('SOX2', 6657),
    ('DMRTA1', 63951),
    ('NDP', 4693),
    ('GNG3', 2785),
    ('NDST3', 9348),
    ('GRIA4', 2893),
    ('SATB2', 23314),
    ('CXCR4', 7852),
    ('CYYR1', 116159),
    ('SCG5', 6447),
    ('ALDH1A3', 220),
    ('KIF26A', 26153),
    ('C6orf117', 112609),
    ('BOC', 91653),
    ('PRLR', 5618),
    ('C4orf18', 51313),
    ('CA14', 23632),
    ('POU3F2', 5454),
    ('SEMA6A', 57556),
    ('IMPG2', 50939),
    ('GABRA5', 2558),
    ('EGFL11', 346007),
    ('NRL', 4901),
    ('MAB21L2', 10586),
    ('TTR', 7276),
    ('NPR3', 4883),
    ('TBR1', 10716),
    ('FSTL5', 56884),
    ('TMEM16B', 57101),
    # ('C5orf23', 4883),  # this is a synonym of NPR3
    ('GNB3', 2784),
    ('DMD', 1756),
    ('PCDH21', 92211),
    ('USH2A', 7399),
    ('RCVRN', 5957),
    ('PDE6H', 5149),
    ('RASGRF2', 5924),
    ('FOXG1B', 2290),
    ('SAMD3', 154075),
    ('TSHZ3', 57616),
    ('MAK', 4117),
    ('PPP2R2B', 5521),
    ('RD3', 343035),
    ('FAM19A4', 151647),
    ('KCNA1', 3736),
    ('EOMES', 8320),
    ('KHDRBS2', 202559),
    ('RBM24', 221662),
    ('UNC5D', 137970),
    ('OAS1', 4938),
    ('GRM8', 2918),
    ('CDH18', 1016),
    ('SNCAIP', 9627),
    ('MPP3', 4356),
    ('CHL1', 10752),
    ('PEX5L', 51555),
    ('LMX1A', 4009),
    ('GPR12', 2835),
    ('FBXL21', 26223),
    ('SH3GL3', 6457),
    ('NID2', 22795),
    ('LINGO2', 158038),
    ('PTHLH', 5744),
    ('CA4', 762),
    ('PRL', 5617),
    ('KCNIP4', 80333),
    ('NEUROD2', 4761),
    ('ST18', 9705),
    ('OTX2', 5015),
])


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
    ('Group D', ALL_NORTHCOTT[ALL_NORTHCOTT.index('KCNA1'):]),
)

NORTHCOTT_GENEID = [(a, [NORTHCOTT_GENEID_MAP[t] for t in b if t in NORTHCOTT_GENEID_MAP]) for a, b in NORTHCOTT_GENES]

NANOSTRING_GENES = (
    ('WNT', ('WIF1', 'TNC', 'GAD1', 'DKK2', 'EMX2'),),
    ('SHH', ('PDLIM3', 'EYA1', 'HHIP', 'ATOH1', 'SFRP1'),),
    ('Group C', ('IMPG2', 'GABRA5', 'EGFL11', 'NRL', 'MAB21L2', 'NPR3'),),  # EYS = EGFL11
    ('Group D', ('KCNA1', 'EOMES', 'KHDRBS2', 'RBM24', 'UNC5D', 'OAS1', 'OTX2')),  # OTX2 added by SB
)

NANOSTRING_GENEID = [(a, [NORTHCOTT_GENEID_MAP[t] for t in b]) for a, b in NANOSTRING_GENES]

SAMPLE_GROUPS_ZHAO = (
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
    ('Healthy cerebellum', (
        'NT1197',
        'NCb1',
        'NCb2',
        'A911105',
        'A508112',
        'A508285',
    ))
)
