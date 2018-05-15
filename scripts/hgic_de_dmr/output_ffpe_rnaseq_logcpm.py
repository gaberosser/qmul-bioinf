from rnaseq import loader, filter, general
from utils import output
import os
import numpy as np


if __name__ == "__main__":
    min_cpm = 1
    obj = loader.load_by_patient('all', type='ffpe')
    samples = [
        'NH15_1661DEF2C',
        'NH15_1877_SP1C',
        'NH15_2101_DEF1A',
        'NH16_270_DEF1Ereplacement',
        'NH16_616DEF1B',
        'NH16_677_SP1A',
        'NH16_1574DEF1A',
        'NH16_1976_DEF1Areplacement',
        'NH16_2063_DEF1Areplacement',
        'NH16_2214DEF1A',
        'NH16_2255DEF1B2',
        'NH16_2806DEF3A1',
    ]

    # remove duplicates
    dat = obj.data.loc[:, samples]
    dat = filter.filter_by_cpm(dat, min_cpm=min_cpm, min_n_samples=1)
    cpm = (dat + 1).divide(dat.sum() + 1, axis=1) * 1e6
    general.add_gene_symbols_to_ensembl_data(cpm)

    outdir = output.unique_output_dir("ffpe_logcpm_values")

    cpm.to_excel(os.path.join(outdir, 'cpm_ffpe_filtered.xlsx'))