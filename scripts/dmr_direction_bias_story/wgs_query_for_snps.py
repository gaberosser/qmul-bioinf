import os
import vcf
import collections
import pandas as pd
from utils import setops, output, log
from settings import DATA_DIR_NON_GIT


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

if __name__ == '__main__':
    """
    Here we are going to iterate over a VCF file containing all variants called in our WGS dataset.
    We are looking specifically for any alterations in or in the vicinity of a pre-defined list of candidate genes
    known to be involved in DNA methylation.

    While we are at it, we may generate some overview plots, too.
    """
    indir = os.path.join(DATA_DIR_NON_GIT, 'wgs', 'x17067/2017-12-12')
    vcf_fn = os.path.join(indir, 'merged.vcf.gz')

    membership = collections.defaultdict(list)
    save_every_n = 1000

    with open(vcf_fn, 'rb') as f:
        reader = vcf.Reader(f, compressed=True)
        for i, rec in enumerate(reader):
            if i % 10000 == 0:
                logger.info("Processed %d VCF entries", i)
            if i % save_every_n == 0:
                for call in rec.samples:
                    if call.called:
                        membership[call.sample].append(str(rec))
