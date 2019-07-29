import os
import gzip
import pickle
import vcf
import vcf.utils
import collections
import pandas as pd
from utils import setops, output, log
from settings import DATA_DIR_NON_GIT
from scripts.hgic_final import consts


logger = log.get_console_logger()


METH_GENES = [
    "A0A096LPK6", "AICDA", "ALKBH1", "ALKBH2", "ALKBH3", "APEX1", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B",
    "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H", "ASZ1", "ATF7IP", "ATRX", "BAZ2A", "BEND3", "BRCA1",
    "CTCF", "CTCFL", "DDX4", "DMAP1", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DPPA3", "EHMT1", "EHMT2", "EZH2",
    "FAM129B", "FKBP6", "FOS", "FTO", "GATA3", "GATAD2A", "GNAS", "GRHL2", "GSK3A", "GSK3B", "H1FOO", "HELLS", "HEMK1",
    "KCNQ1OT1", "KDM1B", "KMT2A", "KMT2E", "MAEL", "MBD1", "MBD2", "MBD3", "MECP2", "METTL4", "MGMT", "MIS18A", "MORC1",
    "MOV10L1", "MPHOSPH8", "MTA2", "MTRR", "MYC", "N6AMT1", "OTUD4", "PARP1", "PICK1", "PIK3CA", "PIKC3A", "PIWIL2",
    "PIWIL4", "PLD6", "PPM1D", "PRDM14", "PRMT5", "PRMT7", "RLF", "SPI1", "STPG4", "TDG", "TDRD1", "TDRD12", "TDRD5",
    "TDRD9", "TDRKH", "TET1", "TET2", "TET3", "TRIM28", "UHRF1", "UHRF2", "USP7", "USP9X", "WT1", "ZFP57", "ZMPSTE24"
]

def run_one_count(res_arr, key1='GIC', key2='iNSC'):
    this_var_count = collections.defaultdict(collections.Counter)
    for t in res_arr:
        if t[key2] is None:
            this_var_count[t[key1].CHROM]['%s only' % key1] += 1
        elif t[key1] is None:
            this_var_count[t[key2].CHROM]['%s only' % key2] += 1
        elif (t[key1].samples[0].gt_nums == '1/1') and (t[key2].samples[0].gt_nums == '0/1'):
            this_var_count[t[key1].CHROM]['%s hom %s het' % (key1, key2)] += 1
        elif (t[key2].samples[0].gt_nums == '1/1') and (t[key1].samples[0].gt_nums == '0/1'):
            this_var_count[t[key1].CHROM]['%s hom %s het' % (key2, key1)] += 1
        elif t[key2].samples[0].gt_nums ==  t[key1].samples[0].gt_nums:
            this_var_count[t[key1].CHROM]['same'] += 1
        else:
            this_var_count[t[key1].CHROM]['other different'] += 1
    return this_var_count

if __name__ == '__main__':
    """
    In the script wgs_query_for_snps we have dumped data to gzipped pickle files containing:
    a) Variants that differ in their call in GIC and iNSC
    b) Variants associated in any way with methylation-related genes
    Now we are going to analyse those outputs.
    """
    pids = consts.PIDS
    contigs = set(['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])

    indir = os.path.join(output.OUTPUT_DIR, 'wgs_query')

    var_count = collections.OrderedDict()
    meth_count = collections.OrderedDict()

    for pid in pids:
        var_fn = os.path.join(indir, "%s_delta_variants.pkl.gz" % pid)
        meth_fn = os.path.join(indir, "%s_meth_variants.pkl.gz" % pid)

        if not os.path.isfile(var_fn):
            raise AttributeError("Unable to find the required input file %s" % var_fn)

        if not os.path.isfile(meth_fn):
            raise AttributeError("Unable to find the required input file %s" % meth_fn)

        # with gzip.open(var_fn, 'rb') as f:
        #     this_var = pickle.load(f)
        #     var_count[pid] = run_one_count(this_var)

        with gzip.open(meth_fn, 'rb') as f:
            this_meth = pickle.load(f)
            meth_count[pid] = run_one_count(this_meth)


