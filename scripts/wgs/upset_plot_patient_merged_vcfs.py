import os
import gzip
from glob import glob
import pickle
import re
import itertools
import vcf.utils
import collections
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils import setops, output, log, dictionary
from plotting import venn
from settings import DATA_DIR
from scripts.hgic_final import consts


logger = log.get_console_logger()


def classify_comparison(rec1, rec2, labels=('1', '2')):
    """
    Given two related variant records, classify (e.g. '1 only').
    Either of the records can be None, indicating absence.
    :param rec1:
    :param rec2:
    :param labels: Optional. These labels are used to
    :return:
    """
    if rec1 is None and rec2 is None:
        return
    elif rec2 is None:
        return "{0} only".format(labels[0])
    elif rec1 is None:
        return "{0} only".format(labels[1])
    elif (rec1.samples[0]['GT'] == '1/1') and (rec2.samples[0]['GT'] == '0/1'):
        return "{0} hom {1} het".format(*labels)
    elif (rec2.samples[0]['GT'] == '1/1') and (rec1.samples[0]['GT'] == '0/1'):
        return "{0} hom {1} het".format(*labels[::-1])
    elif rec1 == rec2:
        return 'same'
    else:
        return 'other'


def process_one_vcf(fn):
    res = {}
    variants_gic_only = []
    variants_hom_het = []
    it = vcf.Reader(filename=fn)

    for i, rec in enumerate(it):
        if i % 50000 == 0:
            logger.info(
                "Processed %d variants.",
                i,
            )
        if i % store_every == 0:
            cl = classify_comparison(gic_rec, insc_rec, labels=('GIC', 'iNSC'))
            if (cl == 'GIC only'):
                variants_gic_only[pid].append(gic_rec)
            elif (cl == 'GIC hom iNSC het'):
                variants_hom_het[pid].append((gic_rec, insc_rec))


if __name__ == '__main__':
    base_indir = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12/')
    meta_fn = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12/', 'sources.csv')
    indir = os.path.join(base_indir, 'paired_vcfs')

    store_every = 20
    pids = consts.PIDS


    variants_hom_het = {}
    variants_gic_only = {}

    for pid in pids:
        fn = os.path.join(indir, "%s.vcf.gz" % pid)
        variants_gic_only[pid] = []
        variants_hom_het[pid] = []
        logger.info("Patient %s", pid)
        it = vcf.Reader(filename=fn)

        for i, (gic_rec, insc_rec) in enumerate(it):
            if i % 50000 == 0:
                logger.info(
                    "Processed %d variants.",
                    i,
                )
            if i % store_every == 0:
                cl = classify_comparison(gic_rec, insc_rec, labels=('GIC', 'iNSC'))
                if (cl == 'GIC only'):
                    variants_gic_only[pid].append(gic_rec)
                elif (cl == 'GIC hom iNSC het'):
                    variants_hom_het[pid].append((gic_rec, insc_rec))

