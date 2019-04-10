import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
import os
import operator
import multiprocessing as mp


import references
from utils import output
from rnaseq import loader


def eval_one_kde(xi):
    # k constant down columns
    # l constant across rows
    k, l = np.meshgrid(xi.values, xi.values)

    # mean over columns = sum over l = sum over kernels
    return stats.poisson.cdf(
        k,
        l + 0.5,
    ).mean(axis=0)


if __name__ == "__main__":
    outdir = output.unique_output_dir()
    obj = loader.load_by_patient(['018', '019', '030', '031'], source='star')
    X = obj.data
    r = 0.5  # offset for Poisson kernels

    # arbitrary gene set: GBM signalling

    gs1 = [
        'ENSG00000077782',
        'ENSG00000133056',
        'ENSG00000136158',
        'ENSG00000110092',
        'ENSG00000169032',
        'ENSG00000139687',
        'ENSG00000171862',
        'ENSG00000140443',
        'ENSG00000142208',
        'ENSG00000105173',
        'ENSG00000147889',
        'ENSG00000140992',
        'ENSG00000169047',
        'ENSG00000150907',
        'ENSG00000116285',
        'ENSG00000145675',
        'ENSG00000165699',
        'ENSG00000065361',
        'ENSG00000177885',
        'ENSG00000011405',
        'ENSG00000139144',
        'ENSG00000105647',
        'ENSG00000121879',
        'ENSG00000051382',
        'ENSG00000171608',
        'ENSG00000105221',
        'ENSG00000117020',
        'ENSG00000197122',
        'ENSG00000109458',
        'ENSG00000105851',
        'ENSG00000118689',
        'ENSG00000184481',
        'ENSG00000134853',
        'ENSG00000113721',
        'ENSG00000141736',
        'ENSG00000146648',
        'ENSG00000066468',
        'ENSG00000105976',
        'ENSG00000110395',
        'ENSG00000196712',
        'ENSG00000213281',
        'ENSG00000133703',
        'ENSG00000174775',
        'ENSG00000078061',
        'ENSG00000157764',
        'ENSG00000132155',
        'ENSG00000124181',
        'ENSG00000067606',
        'ENSG00000065675',
        'ENSG00000027075',
        'ENSG00000126583',
        'ENSG00000163932',
        'ENSG00000163558',
        'ENSG00000166501',
        'ENSG00000154229',
        'ENSG00000197943',
        'ENSG00000126934',
        'ENSG00000034152',
        'ENSG00000065559',
        'ENSG00000137764',
        'ENSG00000108984',
        'ENSG00000076984',
        'ENSG00000100030',
        'ENSG00000102882',
        'ENSG00000111276',
        'ENSG00000123374',
        'ENSG00000124762',
        'ENSG00000103197',
        'ENSG00000123080',
        'ENSG00000147883',
        'ENSG00000110092',
        'ENSG00000118971',
        'ENSG00000118971',
        'ENSG00000135446',
        'ENSG00000105810',
        'ENSG00000101412',
        'ENSG00000147889',
        'ENSG00000135679',
        'ENSG00000198625',
        'ENSG00000141510',
        'ENSG00000100393',
        'ENSG00000149311',
        'ENSG00000012048',
        'ENSG00000139618',
        'ENSG00000012048',
        'ENSG00000116062',
    ]

    # kde for one gene
    counts = np.arange(8000)

    x1 = X.iloc[0]
    g1 = references.ensembl_to_gene_symbol(X.index[0])
    n = float(len(x1))
    p = X.shape[0]
    fr1 = reduce(operator.add, (stats.poisson.pmf(counts, t + r) for t in x1))

    Fr1 = fr1.cumsum() / n

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(counts, fr1, 'k')
    ax2.plot(counts, Fr1, 'r')
    ax1.set_ylabel('Kernels')
    ax2.set_ylabel('KDE')
    ax2.grid('off')
    ax1.set_title("KDE for %s" % g1)

    fig.tight_layout()
    # fig.savefig(os.path.join(outdir, "kde_example_%s.png" % g1), dpi=200)

    # add annotation lines
    for v in x1[:10]:
        ax2.plot([v, v], [0, Fr1[v]], linestyle='--', color='r')
        ax2.plot([v, counts[-1]], [Fr1[v], Fr1[v]], linestyle='--', color='r')

    # fig.savefig(os.path.join(outdir, "kde_example_%s_annotated.png" % g1), dpi=200)

    # run for all genes and samples
    pool = mp.Pool()
    jobs = {}
    for ei, xi in X.iterrows():
        jobs[ei] = pool.apply_async(eval_one_kde, args=(xi,))

    pool.close()
    pool.join()

    Fr = pd.DataFrame(dict([(k, v.get()) for k, v in jobs.items()])).transpose()
    Fr = Fr.loc[X.index]
    Fr.columns = X.columns
