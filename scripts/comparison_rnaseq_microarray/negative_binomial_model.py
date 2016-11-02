from scripts.comparison_rnaseq_microarray import load_illumina_data
from microarray.process import aggregate_by_probe_set
import numpy as np

from scipy.stats import nbinom
from statsmodels.base.model import GenericLikelihoodModel

marray_data, pvals = load_illumina_data.load_normed_microarray_data(pval=None, return_pvals=True)

# reduce to sign genes
marray_data = marray_data.loc[(pvals < 0.05).all(axis=1), :]

probe_set = load_illumina_data.load_illumina_array_library()
marray_ann = load_illumina_data.add_gene_symbol_column(marray_data, probe_set)
marray_by_gene = aggregate_by_probe_set(marray_ann, method="median")

def ll_nbinom(y, X, beta, alph):
    """

    :param y: The responses
    :param X: The regressors
    :param beta: Vector of coefficients
    :param alph: Negative binomial heterogeneity parameter
    :return: Log likelihood
    """
    mu = np.exp(np.dot(X, beta))  # expectation
    size = 1 / float(alph)  # r parameter: number of trials
    prob = size / (size + mu)  # or 1 / (1 + alph * mu): probability of success
    ll = nbinom.logpmf(y, size, prob)
    return ll


class NBin(GenericLikelihoodModel):
     def __init__(self, endog, exog, **kwds):
         super(NBin, self).__init__(endog, exog, **kwds)
     def nloglikeobs(self, params):
         alph = params[-1]
         beta = params[:-1]
         ll = ll_nbinom(self.endog, self.exog, beta, alph)
         return -ll
     def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
         if start_params == None:
             # Reasonable starting values
             start_params = np.append(np.zeros(self.exog.shape[1]), .5)
             start_params[0] = np.log(self.endog.mean())
         return super(NBin, self).fit(start_params=start_params,
                                      maxiter=maxiter, maxfun=maxfun,
                                      **kwds)

M, N = marray_by_gene.shape  # no. genes, no. patients
y = marray_by_gene.values.flatten('C')  # length N X M, flattened by ROWS
X = np.zeros((M * N, M), dtype=np.uint8)

for i in range(M):
    X[(N * i):((N + 1) * i), i] = 1.

mod = NBin(y, X)