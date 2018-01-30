from load_data import rnaseq_data
import pandas as pd
import os
from rnaseq import general
from utils import output
import references


if __name__ == "__main__":
    gois = [
        'BMI1'
    ]
    pids = 'all'
    cell_types = ['GBM']
    unit = 'cpm'
    type = 'cell_culture'
    tax_id = 9606

    if isinstance(cell_types, str):
        cell_types = [cell_types]

    if unit == 'tpm':
        dat = rnaseq_data.load_salmon_by_patient_id(pids, include_control=False, type=type)
        dat = general.ensembl_transcript_quant_to_gene(dat, tax_id=tax_id)
    elif unit == 'cpm':
        obj = rnaseq_data.load_by_patient(pids, type=type, source='star', annotate_by='Ensembl Gene ID', include_control=False)
        dat = obj.data
        dat = dat.divide(dat.sum(axis=0), axis=1) * 1e6
    else:
        raise NotImplementedError()

    if cell_types is not None:
        idx = reduce(
            lambda x, y: x | y,
            [dat.columns.str.contains(t) for t in cell_types],
        )
        dat = dat.loc[:, idx]

    lookup = references.gene_symbol_to_ensembl(gois, tax_id=tax_id)

    vals = dat.loc[lookup]