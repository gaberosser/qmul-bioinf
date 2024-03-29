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

METHYLATION_DIRECTION_COLOURS = {
    # 'hypo': '#89CD61',
    'hypo': '#c70039',
    # 'hyper': '#FF381F',
    'hyper': '#3d3d6b',
}

PATIENT_COLOURS = {
    '018': '#ccffcc',
    '019': '#4dff4d',
    '030': '#00cc00',
    '031': '#004d00',
    '017': '#ffcccc',
    '050': '#ff4d4d',
    '054': '#cc0000',
    '061': '#660000',
    '026': '#ff80ff',
    '052': '#800080'
}

FFPE_RNASEQ_SAMPLES = [
    'NH15_1661DEF2C',
    'NH15_1877_SP1C',
    'NH15_2101_DEF1A',
    'NH16_270_DEF1Ereplacement',
    'NH16_616DEF1B',
    'NH16_677_SP1A',
    'NH16_2063_DEF1Areplacement',
    'NH16_2214DEF1A',
    'NH16_2255DEF1B2',
    'NH16_2806DEF3A1'
]

FFPE_RNASEQ_SAMPLES_ALL = FFPE_RNASEQ_SAMPLES + [
    'NH16_1574DEF1A',
    'NH16_1976_DEF1Areplacement',
]

S1_RNASEQ_SAMPLES_GIC = [
    'GBM018_P12',
    'GBM018_P10',
    'GBM019_P4',
    'GBM019_P3n6',
    'GBM030_P9n10',
    'GBM030_P5',
    'GBM031_P7',
    'GBM031_P4',
    'GBM017_P3',
    'GBM017_P4',
    'GBM050_P7n8',
    'GBM050_P9',
    'GBM054_P4',
    'GBM054_P6',
    'GBM061_P3',
    'GBM061_P5',
    'GBM026_P8',
    'GBM026_P3n4',
    'GBM052_P6n7',
    'GBM052_P4n5',
]

S1_RNASEQ_SAMPLES_INSC = [
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'DURA019_NSC_N8C_P2',
    'DURA019_NSC_N5C1_P2',
    'DURA030_NSC_N16B6_P1',
    'DURA030_NSC_N9_P2',
    'DURA031_NSC_N44B_P2',
    'DURA031_NSC_N44_P3',
    'DURA017_NSC_N3C5_P4',
    'DURA050_NSC_N12_P3',
    'DURA050_NSC_N16_P4',
    'DURA054_NSC_N3C_P2',
    'DURA054_NSC_N2E_P1',
    'DURA061_NSC_N4_P2',
    'DURA061_NSC_N1_P3',
    'DURA026_NSC_N31D_P5',
    'DURA052_NSC_N4_P3',
    'DURA052_NSC_N5_P2',
]

S1_RNASEQ_SAMPLES_FB = [
    'DURA018_FB_P6',
    'DURA019_FB_P7',
    'DURA030_FB_P8',
    'DURA031_FB_P7',
    'DURA017_FB_P7',
    'DURA050_FB_P7',
    'DURA054_FB_P5',
    'DURA061_FB_P7',
    'DURA026_FB_P8',
    'DURA052_FB_P6'
]

S1_RNASEQ_SAMPLES_IPSC = [
    'DURA019_IPSC_N8_P13',
    'DURA030_IPSC_N16B6_P13',
    'DURA031_IPSC_N44B_P10',
    'DURA050_IPSC_N12_P5',
    'DURA054_IPSC_N3C_P11',
    'DURA061_IPSC_N4_P5'
]

S1_RNASEQ_SAMPLES_IAPC = [
    'DURA019_ASTRO_N8C_D20_B2',
    'DURA019_ASTRO_N8C_D20_B3',
    'DURA031_ASTRO_N44B_D20_B1',
    'DURA031_ASTRO_N44B_D20_B2',
    'DURA050_ASTRO_N12_D20_B1',
    'DURA050_ASTRO_N12_D20_B2',
    'DURA052_ASTRO_N8_D20_B2',
    'DURA052_ASTRO_N8_D20_B3'
]

# S1_RNASEQ_SAMPLES_IOPC = []

S1_RNASEQ_SAMPLES = S1_RNASEQ_SAMPLES_GIC + S1_RNASEQ_SAMPLES_INSC

ALL_RNASEQ_SAMPLES = S1_RNASEQ_SAMPLES + [
    'GIBCO_NSC_P4',
    'H9_NSC_1',
    'H9_NSC_2',
]

S1_METHYL_SAMPLES_GIC = [
    'GBM018_P12',
    'GBM018_P10',
    'GBM019_P4',
    'GBM019_P3n6',
    'GBM030_P9',
    'GBM030_P5',
    'GBM031_P7',
    'GBM031_P4',
    'GBM017_P3',
    'GBM017_P4',
    'GBM050_P7n8',
    'GBM050_P9',
    'GBM054_P4',
    'GBM054_P6',
    'GBM061_P3',
    'GBM061_P5',
    'GBM026_P8',
    'GBM026_P3n4',
    'GBM052_P6n7',
    'GBM052_P4n5',
]

# NB: JB has also sequenced some replicates of these iNSC lines, but we don't currently use them
S1_METHYL_SAMPLES_INSC = [
    'DURA018_NSC_N4_P4',
    'DURA018_NSC_N2_P6',
    'DURA019_NSC_N8C_P2',
    'DURA019_NSC_N5C1_P2',
    'DURA030_NSC_N16B6_P1',
    'DURA030_NSC_N9_P2',
    'DURA031_NSC_N44B_P2',
    'DURA031_NSC_N44F_P3',
    'DURA017_NSC_N3C5_P4',
    'DURA050_NSC_N12_P3',
    'DURA050_NSC_N16_P4',
    'DURA054_NSC_N3C_P2',
    'DURA054_NSC_N2E_P1',
    'DURA061_NSC_N4_P2',
    'DURA061_NSC_N1_P3n4',
    'DURA026_NSC_N31D_P5',
    'DURA052_NSC_N4_P3',
    'DURA052_NSC_N5_P2',
]

S1_METHYL_SAMPLES = S1_METHYL_SAMPLES_GIC + S1_METHYL_SAMPLES_INSC

S1_METHYL_SAMPLES_FB = [
    'DURA019_FB_P7',
    'DURA030_FB_P8',
    'DURA031_FB_P7',
    'DURA017_FB_P7',
    'DURA050_FB_P7',
    'DURA054_FB_P5',
    'DURA018_NH15_1877_P6_15/05/2017',
    'DURA026_NH16_270_P8_15/05/2017',
    'DURA052_NH16_2214_P6_14/04/2017'
]

S1_METHYL_SAMPLES_IOPC = [
    'DURA019_IOPC_N8C_B2',
    'DURA031_IOPC_N44B_B2',
    'DURA050_IOPC_N12_B2',
    'DURA052_IOPC_N4_B2',
    'DURA031_OPC',
    'DURA019_OPC',
    'DURA050_OPC',
    'DURA052_OPC',
]

S1_METHYL_SAMPLES_IAPC = [
    'DURA019_IAPC_D20_N8C_B2',
    'DURA019_IAPC_D20_N8C_B3',
    'DURA031_IAPC_D20_N44B_B1',
    'DURA031_IAPC_D20_N44B_B2',
    'DURA050_IAPC_D20_N12_B1',
    'DURA050_IAPC_D20_N12_B2',
    'DURA052_IAPC_D20_N4_B2',
    'DURA052_IAPC_D20_N4_B3',
]

S1_METHYL_SAMPLES_IPSC = [
    'DURA019_IPSC_N8C_P13',
    'DURA030_IPSC_N16B6_P13',
    'DURA031_IPSC_N44B_P10',
    'DURA050_IPSC_N12_P5',
    'DURA054_IPSC_N3C_P11',
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