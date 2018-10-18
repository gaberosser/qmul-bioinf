from rnaseq import gsea
import pandas as pd
from . import consts
from utils.log import get_console_logger
import logging
logger = get_console_logger("signature_classifier")


class ssGSEAClassifier(object):
    """
    Basic classifier that uses pre-defined signatures to score samples and assess classification.
    """
    def __init__(self, signature_dict, **ssgsea_kwds):
        """
        
        :param signature_dict: Dictionary. Keys are the class name, values are iterables of genes / probes
        or any other row index
        :param ssgsea_kwds: Any additional kwargs are passed directly to the ssgsea algorithm.
        """
        self.signatures = signature_dict
        self.ssgsea_kwds = ssgsea_kwds
        # it's useful to maintain a copy of all the signature IDs for validation purposes
        self.all_ids = reduce(lambda x, y: x.union(y), self.signatures.values(), set())

    def score(self, sample_data, check_overlap=True):
        """
        
        :param sample_data: Pandas Series (single sample) or DataFrame (multiple samples) to be classified.
        :param check_overlap: If True (default), check whether all the signature IDs are in the data. If
        any are missing, raise an error.
        :return: 
        """
        # if a single sample has been supplied, create a DataFrame with a single column
        series_in = False
        if isinstance(sample_data, pd.Series):
            sample_data = pd.DataFrame(sample_data)
            series_in = True

        # check that the index is unique
        if sample_data.index.duplicated().any():
            if sample_data.index[sample_data.index.duplicated()].isin(self.all_ids).any():
                raise ValueError("Data index contains duplicate entries that intersect the signature(s).")
            else:
                logger.warn("Data index contains duplicate entries, but these are not in the signatures.")

        # check that the data index contains the signature genes
        if self.all_ids.intersection(sample_data.index) != self.all_ids:
            missing = self.all_ids.difference(sample_data.index)
            if check_overlap:
                raise ValueError("Some signature IDs are not in the data: %s" % ', '.join(missing))
            else:
                logger.warn("Some signature IDs are not in the data: %s", ', '.join(missing))

        # compute score
        ## TODO: check this - syntax for ssgsea has changed
        scores = pd.DataFrame(index=self.signatures.keys(), columns=sample_data.columns)
        for col, dat in sample_data.iteritems():
            scores.loc[:, col] = pd.Series(gsea.ssgsea(dat, self.signatures, **self.ssgsea_kwds))
        # for s, arr in self.signatures.items():
        #     for col, dat in sample_data.iteritems():
        #         scores.loc[s, col] = gsea.ssgsea(dat, arr, **self.ssgsea_kwds)

        if series_in:
            return scores.iloc[:, 0]
        else:
            return scores


class WangGBMClassifier(ssGSEAClassifier):
    def __init__(self, annotate_by='Ensembl Gene ID', **kwargs):
        if annotate_by == 'Ensembl Gene ID':
            super(WangGBMClassifier, self).__init__(consts.WANG_SIGNATURES_BY_ENSEMBL, **kwargs)
        elif annotate_by == 'Approved Symbol':
            super(WangGBMClassifier, self).__init__(consts.WANG_SIGNATURES_BY_SYMBOL, **kwargs)
        else:
            raise AttributeError("Unrecognised annotation type %s." % annotate_by)


"""
TODO now:
Load TCGA data 
- microarray 
- RNASeq
Define the 
"""