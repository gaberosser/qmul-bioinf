SUBGROUPS = {
    'RTK I': ['018', '019', '030', '031'],
    'RTK II': ['017', '050', '054', '061'],
    'MES': ['026', '052']
}
PIDS = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052']

# for some scripts, we may wish to include more patients (e.g. when characterising iPSC lines)
ALL_PIDS = ['018', '019', '030', '031', '017', '050', '054', '061', '026', '052', '044', '049']

SUBGROUP_SET_COLOURS = {
    'RTK I full': '#0d680f',
    'RTK II full': '#820505',
    'MES full': '#7900ad',
    'RTK I partial': '#6ecc70',
    'RTK II partial': '#d67373',
    'MES partial': '#cc88ea',
    'Expanded core': '#4C72B0',
    'Specific': '#f4e842',
    }

S1_RNASEQ_SAMPLES = [
    'GBM018_P12',
    'GBM018_P10',
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'GBM019_P4',
    'GBM019_P3n6',
    'DURA019_NSC_N8C_P2',
    'DURA019_NSC_N5C1_P2',
    'GBM030_P9n10',
    'GBM030_P5',
    'DURA030_NSC_N16B6_P1',
    'DURA030_NSC_N9_P2',
    'GBM031_P7',
    'GBM031_P4',
    'DURA031_NSC_N44B_P2',
    'DURA031_NSC_N44_P3',
    'GBM017_P3',
    'GBM017_P4',
    'DURA017_NSC_N3C5_P4',
    'GBM050_P7n8',
    'GBM050_P9',
    'DURA050_NSC_N12_P3',
    'DURA050_NSC_N16_P4',
    'GBM054_P4',
    'GBM054_P6',
    'DURA054_NSC_N3C_P2',
    'DURA054_NSC_N2E_P1',
    'GBM061_P3',
    'GBM061_P5',
    'DURA061_NSC_N4_P2',
    'DURA061_NSC_N1_P3',
    'GBM026_P8',
    'GBM026_P3n4',
    'DURA026_NSC_N31D_P5',
    'GBM052_P6n7',
    'GBM052_P4n5',
    'DURA052_NSC_N4_P3',
    'DURA052_NSC_N5_P2',
]
ALL_RNASEQ_SAMPLES = S1_RNASEQ_SAMPLES + [
    'GIBCO_NSC_P4',
    'H9_NSC_1',
    'H9_NSC_2',
]

S1_METHYL_SAMPLES = [
    'GBM018_P12',
    'GBM018_P10',
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'GBM019_P4',
    'GBM019_P3n6',
    'DURA019_NSC_N8C_P2',
    'DURA019_NSC_N5C1_P2',
    'GBM030_P9',
    'GBM030_P5',
    'DURA030_NSC_N16B6_P1',
    'DURA030_NSC_N9_P2',
    'GBM031_P7',
    'GBM031_P4',
    'DURA031_NSC_N44B_P2',
    'DURA031_NSC_N44F_P3',
    'GBM017_P3',
    'GBM017_P4',
    'DURA017_NSC_N3C5_P4',
    'GBM050_P7n8',
    'GBM050_P9',
    'DURA050_NSC_N12_P3',
    'DURA050_NSC_N16_P4',
    'GBM054_P4',
    'GBM054_P6',
    'DURA054_NSC_N3C_P2',
    'DURA054_NSC_N2E_P1',
    'GBM061_P3',
    'GBM061_P5',
    'DURA061_NSC_N4_P2',
    'DURA061_NSC_N1_P3n4',
    'GBM026_P8',
    'GBM026_P3n4',
    'DURA026_NSC_N31D_P5',
    'GBM052_P6n7',
    'GBM052_P4n5',
    'DURA052_NSC_N4_P3',
    'DURA052_NSC_N5_P2',
]

ALL_METHYL_SAMPLES = S1_METHYL_SAMPLES + [
    'GIBCONSC_P4',
    'H9 NPC 1',
    'H9 NPC 2',
]

# parameters

DE_PARAMS = {
    'lfc': 1,
    'fdr': 0.01,
    'method': 'QLGLM'
}

DMR_PARAMS = {
    'd_max': 400,
    'n_min': 6,
    'delta_m_min': 1.4,
    'alpha': 0.01,
    'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
    'test_kwargs': {},
    'n_jobs': None  # NB: fill this in within the script
}