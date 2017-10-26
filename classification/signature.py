from rnaseq import gsea
import pandas as pd
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
        scores = pd.DataFrame(index=self.signatures.keys(), columns=sample_data.columns)
        for s, arr in self.signatures.items():
            scores.loc[s] = gsea.ssgsea(sample_data, arr, **self.ssgsea_kwds)

        if series_in:
            return scores.iloc[:, 0]
        else:
            return scores


WANG_SIGNATURES_BY_SYMBOL = {
    'CL': [
        "HN1",
        "RAB33A",
        "HDAC2",
        "MYT1",
        "MTSS1",
        "HOXD3",
        "GPR17",
        "PTTG1",
        "KLRC3",
        "HRASLS",
        "TCP1",
        "NPPA",
        "PFDN2",
        "CA10",
        "EPHB1",
        "UGT8",
        "PAK7",
        "SLC1A1",
        "NARF",
        "DCTN3",
        "SMPD3",
        "ZNF804A",
        "RASL11B",
        "MYB",
        "PDGFRA",
        "ERBB3",
        "CLGN",
        "SOX10",
        "BCL11A",
        "NMU",
        "ZNF643",
        "CDKN1C",
        "JPH3",
        "PCDHA9",
        "IL1RAPL1",
        "MAST1",
        "VIPR2",
        "SIM2",
        "BAMBI",
        "PKMYT1",
        "PLCB4",
        "SLC17A6",
        "KLRK1",
        "CENPJ",
        "NHLH1",
        "GABRB3",
        "KLRC4",
        "KCNK3",
        "GRID2",
        "DACH1",
    ],
    'PN': [
        "PTPRA",
        "ELOVL2",
        "MLC1",
        "SOX9",
        "ARNTL",
        "DENND2A",
        "BBS1",
        "ABLIM1",
        "PAX6",
        "ZHX3",
        "USP8",
        "PLCG1",
        "CDH4",
        "RASGRP1",
        "ACSBG1",
        "CST3",
        "BCKDHB",
        "LHFP",
        "VAV3",
        "ACSL3",
        "EYA2",
        "SEPT11",
        "SLC4A4",
        "SLC20A2",
        "C14orf159",
        "CTNND1",
        "ZFHX4",
        "SPRY2",
        "ZNF45",
        "NCOA1",
        "PLCE1",
        "DTNA",
        "POLRMT",
        "SALL1",
        "TYK2",
        "TJP1",
        "MEOX2",
        "FGFR3",
        "STXBP3",
        "GRIK1",
        "GATM",
        "UPF1",
        "NPEPL1",
        "KIAA0494",
        "RBCK1",
        "PHKB",
        "SLC3A2",
        "PPARGC1A",
        "PNPLA6",
        "MYO5C",
    ],
    'MES': [
        "ARPC1B",
        "S100A11",
        "CTSC",
        "GLIPR1",
        "NNMT",
        "VDR",
        "RGS2",
        "CTSB",
        "TGFBI",
        "PLAUR",
        "LY96",
        "BCL3",
        "TNFAIP8",
        "IER3",
        "PRSS23",
        "IL7R",
        "RAB27A",
        "RUNX1",
        "P4HA2",
        "CYP1B1",
        "BACE2",
        "ACPP",
        "FTL",
        "SLPI",
        "RAC2",
        "RARRES1",
        "SYNGR2",
        "THBS1",
        "IL6",
        "CAV1",
        "PI3",
        "CDCP1",
        "ITGB1",
        "LOX",
        "CD72",
        "COL1A2",
        "ANPEP",
        "MMP7",
        "SPAG4",
        "BNC2",
        "NDRG1",
        "CNN2",
        "LUM",
        "PTGS2",
        "COL3A1",
        "COL5A1",
        "SDC1",
        "COL1A1",
        "GPRC5A",
        "COL15A1",
    ],
}


class WangGBMClassifier(ssGSEAClassifier):
    def __init__(self, **kwargs):
        super(WangGBMClassifier, self).__init__(WANG_SIGNATURES_BY_SYMBOL, **kwargs)


"""
TODO now:
Load TCGA data 
- microarray 
- RNASeq
Define the 
"""