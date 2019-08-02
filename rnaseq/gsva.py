import multiprocessing as mp
import operator
from functools import partial

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

from rnaseq import loader
from utils import output, log, rinterface, reference_genomes

logger = log.get_console_logger()


if rinterface.RFUNCTIONS_PRESENT and rinterface.RPANDAS_PRESENT:
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri, r
    pandas2ri.activate()

    def _run_gsea(df, genesets, method='ssgsea', verbose=False, **kwargs):
        rdata = r('as.matrix')(df)
        rgenesets = robjects.ListVector(genesets)
        res = r('gsva')(rdata, rgenesets, method=method, verbose=verbose, **kwargs)
        py_res = pandas2ri.ri2py_dataframe(res)
        py_res.index = r('rownames')(res)
        # py_res.columns = r('colnames')(res)
        py_res.columns = df.columns
        return py_res

    run_gsea = rinterface.RFunctionDeferred(_run_gsea, imports=['GSVA'])
    ssgsea = rinterface.RFunctionDeferred(partial(_run_gsea, method='ssgsea'), imports=['GSVA'])
    gsva = rinterface.RFunctionDeferred(partial(_run_gsea, method='gsva'), imports=['GSVA'])

def eval_one_kde_poisson(xi):
    # k constant down columns
    # l constant across rows
    k, l = np.meshgrid(xi.values, xi.values)

    # mean over columns = sum over l = sum over kernels
    return stats.poisson.cdf(
        k,
        l + 0.5,
    ).mean(axis=0)


def eval_one_kde_gaussian(xi):
    return xi.rank(method='average') / float(xi.size)


def compute_ks_statistic(z_j, g_k, tau=1.0):
    p = len(z_j)
    r_j = np.abs(p * .5 - z_j)

    ranked_ix = z_j.sort_values().index
    g_in = ranked_ix.isin(g_k)
    r_j_ranked = r_j.loc[ranked_ix]

    a = (r_j_ranked ** tau * g_in).cumsum()
    b = a[-1]
    c = (~g_in).cumsum().astype(float)
    d = p - g_in.sum()

    return a/b, c/d


def compute_ks_statistic_no_norm(z_j, g_k, tau=1.0):
    ranked_ix = z_j.sort_values().index
    g_in = ranked_ix.isin(g_k)
    z_j_ranked = z_j.loc[ranked_ix]

    a = (z_j_ranked ** tau * g_in).cumsum()
    b = a[-1]
    c = (~g_in).cumsum().astype(float)
    d = p - g_in.sum()

    return a/b, c/d


class GSVA(object):
    def __init__(self, X, gene_sets=None, tau=1., njob=None):
        self.X = X
        self.tau = tau
        self.z_ij = None
        self.__gene_sets = {}
        if gene_sets is not None:
            self.gene_sets = gene_sets
        self.run_normalisation(njob=njob)

    @property
    def n(self):
        return self.X.shape[1]

    @property
    def p(self):
        return self.X.shape[0]

    @property
    def gene_sets(self):
        return self.__gene_sets

    @gene_sets.setter
    def gene_sets(self, gs_dict):
        self.__gene_sets = {}
        for k, v in gs_dict.items():
            the_ix = pd.Index(v).intersection(self.X.index)
            if len(the_ix) < len(v):
                logger.warn(
                    "%d genes in gene set %s (original size: %d) were not found in the count matrix and dropped.",
                    len(v) - len(the_ix),
                    k,
                    len(v)
                )
            if the_ix.duplicated().any():
                logger.warn(
                    "Gene set %s contains %d duplicates, which will be dropped.",
                    k,
                    the_ix.duplicated().sum()
                )
                the_ix = the_ix[~the_ix.duplicated()]
            self.__gene_sets[k] = the_ix

    @property
    def kde_func(self):
        raise NotImplementedError("KDE implementation is included in the derived classes - please run one of those.")

    def run_normalisation(self, njob=None):
        if njob is None:
            njob = mp.cpu_count()

        if njob > 1:
            pool = mp.Pool()
            jobs = {}
            for ei, xi in self.X.iterrows():
                jobs[ei] = pool.apply_async(self.kde_func, args=(xi,))

            pool.close()
            pool.join()

            Fr = pd.DataFrame(dict([(k, v.get()) for k, v in jobs.items()])).transpose()
            Fr = Fr.loc[self.X.index]
            Fr.columns = self.X.columns

        else:
            Fr = pd.DataFrame(dict([
                                       (k, eval_one_kde_poisson(xi)) for k, xi in self.X.iterrows()
                                       ]))

        self.Fr = Fr
        # convert to normalized ranked data (genes ranked within each sample)
        self.z_ij = Fr.rank(axis=0, method='average')

    def plot_enrichment_one_gene_set_one_sample(self, sample_name=None, gene_set=None, add_scores=True):
        if sample_name is None:
            sample_name = self.X.columns[0]

        if gene_set is None:
            gene_set = self.gene_sets.keys()[0]

        if gene_set not in self.gene_sets:
            raise KeyError("gene_set must be in self.gene_sets")

        if sample_name not in self.X.columns:
            raise KeyError("sample_name must be in the columns of self.X")

        a, b = compute_ks_statistic(self.z_ij[sample_name], self.gene_sets[gene_set], tau=self.tau)

        fig, ax = plt.subplots()
        ax.plot(range(1, self.p + 1), a, c='k', label='Observed for %s' % sample_name)
        ax.plot(range(1, self.p + 1), b, c='b', linestyle='--', label='Null')
        if add_scores:
            # add ES1 line
            pos = np.abs(a - b).values.argmax()
            y0 = min(a[pos], b[pos])
            y1 = max(a[pos], b[pos])
            ax.plot([pos + 1, pos + 1], [y0, y1], c='r', linestyle=':', zorder=10, label='ES1')
            # add ES2 lines
            pos_max = (a - b).values.argmax()
            y_max_0 = min(a[pos_max], b[pos_max])
            y_max_1 = max(a[pos_max], b[pos_max])
            pos_min = (a - b).values.argmin()
            y_min_0 = min(a[pos_min], b[pos_min])
            y_min_1 = max(a[pos_min], b[pos_min])
            ax.plot([pos_min + 1, pos_min + 1], [y_min_0, y_min_1], c='y', label='ES2(-)')
            ax.plot([pos_max + 1, pos_max + 1], [y_max_0, y_max_1], c='g', label='ES2(+)')

        ax.set_xlabel('Rank (ascending)')
        ax.set_ylabel('KS statistic')
        ax.legend()
        fig.tight_layout()

        return fig, ax

    def run_enrichment(self, gene_sets=None):
        """
        Run the final step of the GSVA algorithm, in which enrichment scores are computed.
        We compute three variants, all of which are discussed briefly in the paper and user guide.
        :param gene_sets: If supplied, this is an iterable containing the names of gene sets to test. Otherwise, they
        are all tested.
        :return: Dictionary with the following entries
        v_jk: The KS random walk statistic, from which all enrichment scores are computed.
        es1_jk: ES1 is the traditional KS statistic, defined as the maximum absolute value of v_ij. This works for
        gene signatures that go up or down.
        es2_jk: ES2 is the alternative ES described in the paper, taken as the difference between the maximum positive
        deviation and the (absolute value of the) maximum negative deviation. This is designed for signatures with
        genes that all go in the same direction; it penalises cases where they appear in both tails.
        es3_jk: As for ES2, but we sum the two values, which is identical to the Kuiper test. This considers appearance
        in either tail (or a mixture) as an enrichment.
        """
        ## TODO: might consider making this parallel for larger jobs
        if len(self.gene_sets) == 0:
            raise AttributeError("Must set gene_sets before running enrichment analysis.")

        if self.z_ij is None:
            raise AttributeError("Method run_normalisation() has not been called. Do so before proceeding.")

        es1_jk = {}
        es2_jk = {}
        es3_jk = {}
        v_jk = {}
        for k, gk in self.gene_sets.iteritems():
            if gene_sets is not None and k not in gene_sets:
                continue

            v_j = {}
            for col in self.X.columns:
                a, b = compute_ks_statistic(self.z_ij[col], gk, tau=self.tau)
                v_j[col] = a - b

            v_j = pd.DataFrame(v_j)
            v_jk[k] = v_j

            # ES (1): maximum deviation from zero
            # This has a bimodal distribution under the null and must be tested from both sides (?)
            # use this if the gene set contains both up- and downregulated genes
            es1_jk[k] = np.abs(v_j).max(axis=0)

            # ES (2): Alternative with Gaussian distribution under the null
            a = v_j.max(axis=0)
            a[a < 0] = 0.
            b = v_j.min(axis=0)
            b[b > 0] = 0.

            es2_jk[k] = a + b  # b is negative
            es3_jk[k] = a - b  # b is negative

        v_jk = v_jk
        es1_jk = pd.DataFrame(es1_jk)
        es2_jk = pd.DataFrame(es2_jk)
        es3_jk = pd.DataFrame(es3_jk)

        return {
            'v_jk': v_jk,
            'es1_jk': es1_jk,
            'es2_jk': es2_jk,
            'es3_jk': es3_jk,
        }


class GSVAForCounts(GSVA):
    @property
    def kde_func(self):
        return eval_one_kde_poisson

    def run_normalisation(self, njob=None):
        if njob is None:
            njob = mp.cpu_count()

        if njob > 1:
            pool = mp.Pool()
            jobs = {}
            for ei, xi in self.X.iterrows():
                jobs[ei] = pool.apply_async(self.kde_func, args=(xi,))

            pool.close()
            pool.join()

            Fr = pd.DataFrame(dict([(k, v.get()) for k, v in jobs.items()])).transpose()
            Fr = Fr.loc[self.X.index]
            Fr.columns = self.X.columns

        else:
            Fr = pd.DataFrame(dict([
                                       (k, eval_one_kde_poisson(xi)) for k, xi in self.X.iterrows()
                                       ]))

        self.Fr = Fr
        # convert to normalized ranked data (genes ranked within each sample)
        self.z_ij = Fr.rank(axis=0, method='average')

    def plot_kde_one_gene(self, gene, n_annot=10, gene_ttl=None):
        """
        Generate a plot showing the Poisson KDE for the specified gene, optionally adding lookup lines
        :param gene: The gene of interest. Must be in the index of self.X
        :param n_annot: Add this many lines (chosen as the first n_annot genes) to the plot to show the norming process.
        To disable annotations, set this to 0.
        :param gene_ttl: If supplied, this is an alternative title for the gene (e.g. gene symbol rather than ID)
        :return:
        """
        if gene_ttl is None:
            gene_ttl = gene
        x = self.X.loc[gene]
        counts = np.arange(int(x.max() * 1.1))
        fr1 = reduce(operator.add, (stats.poisson.pmf(counts, t + 0.5) for t in x))
        Fr1 = fr1.cumsum() / float(self.n)

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        ax1.plot(counts, fr1, 'k')
        ax2.plot(counts, Fr1, 'r')
        ax1.set_ylabel('Kernels')
        ax2.set_ylabel('KDE')
        ax2.grid('off')
        ax1.set_title("KDE for %s" % gene_ttl)

        fig.tight_layout()

        if n_annot > 0:
            # add annotation lines
            for v in x[:10]:
                ax2.plot([v, v], [0, Fr1[v]], linestyle='--', color='r')
                ax2.plot([v, counts[-1]], [Fr1[v], Fr1[v]], linestyle='--', color='r')

        return fig, ax1, ax2


class GSVAForNormedData(GSVA):
    """
    Run GSVA on gene expression data that have been normalised already, such a FPKM or TPM.
    NB for CPM, it might be better to use the counts variant as the KDE is a little more appropriate?
    """
    @property
    def kde_func(self):
        return eval_one_kde_gaussian

    def run_normalisation(self, njob=None):
        self.Fr = self.X.rank(axis=1, method='average') / float(self.n)
        self.z_ij = self.Fr.rank(axis=0, method='average')

    def plot_kde_one_gene(self, gene, n_annot=10, gene_ttl=None):
        raise NotImplementedError("TODO: refactor this part")



if __name__ == "__main__":
    outdir = output.unique_output_dir()
    obj = loader.load_by_patient(['018', '019', '030', '031'], source='star')
    X = obj.data
    r = 0.5  # offset for Poisson kernels
    tau = 1.  # weighting in KS statistic

    # arbitrary gene set: GBM signalling

    gk = [
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
        'ENSG00000118971',
        'ENSG00000135446',
        'ENSG00000105810',
        'ENSG00000101412',
        'ENSG00000135679',
        'ENSG00000198625',
        'ENSG00000141510',
        'ENSG00000100393',
        'ENSG00000149311',
        'ENSG00000012048',
        'ENSG00000139618',
        'ENSG00000116062',
    ]

    # kde for one gene
    counts = np.arange(8000)

    example_col = X.columns[0]
    example_ens = X.index[0]
    example_gene = reference_genomes.ensembl_to_gene_symbol(example_ens)

    x1 = X.loc[example_ens]
    n = float(len(x1))
    p = X.shape[0]
    fr1 = reduce(operator.add, (stats.poisson.pmf(counts, t + r) for t in x1))

    Fr1 = fr1.cumsum() / n

    # run for all genes (rows)
    pool = mp.Pool()
    jobs = {}
    for ei, xi in X.iterrows():
        jobs[ei] = pool.apply_async(eval_one_kde_poisson, args=(xi,))

    pool.close()
    pool.join()

    Fr = pd.DataFrame(dict([(k, v.get()) for k, v in jobs.items()])).transpose()
    Fr = Fr.loc[X.index]
    Fr.columns = X.columns

    # convert to normalized ranked data (genes ranked within each sample)
    z_ij = Fr.rank(axis=0, method='average')
    r_ij = np.abs(p * .5 - z_ij)

    v_j = {}
    for col in X.columns:
        a, b = compute_ks_statistic(z_ij[col], gk, tau=tau)
        v_j[col] = a - b

    v_j = pd.DataFrame(v_j)

    # ES (1): maximum deviation from zero
    # This has a bimodal distribution under the null and must be tested from both sides (?)
    # use this if the gene set contains both up- and downregulated genes
    es1_j = np.abs(v_j).max(axis=0)

    # ES (2): Alternative with Gaussian distribution under the null
    a = v_j.max(axis=0)
    a[a < 0] = 0.
    b = v_j.min(axis=0)
    b[b > 0] = 0.

    es2_j = a + b