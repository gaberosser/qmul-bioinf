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

if __name__ == '__main__':
    """
    Here we are going to iterate over a VCF file containing all variants called in our WGS dataset.
    We are looking specifically for any alterations in or in the vicinity of a pre-defined list of candidate genes
    known to be involved in DNA methylation.

    While we are at it, we may generate some overview plots, too.
    """
    pids = consts.PIDS
    contigs = set(['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])

    base_indir = os.path.join(DATA_DIR_NON_GIT, 'wgs', 'x17067/2017-12-12')
    meta_fn = os.path.join(base_indir, 'sources.csv')

    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    meta.loc[:, 'patient_id'] = ["%03d" % t for t in meta.patient_id]

    outdir = output.unique_output_dir()

    for pid in pids:
        logger.info("Patient %s", pid)
        this_meta = meta.loc[meta.patient_id == pid]
        in_fns = []
        readers = {}
        for t in this_meta.index:
            the_fn = os.path.join(base_indir, t, "%s.vcf.gz" % meta.loc[t, 'sample'])
            in_fns.append(the_fn)
            readers[this_meta.loc[t, 'type']] = vcf.Reader(filename=the_fn, compressed=True)
        logger.info("Found %d conditions to compare: %s", len(readers), ', '.join(readers.keys()))

        it = vcf.utils.walk_together(*readers.values())

        this_res = []
        this_res_meth = []

        for i, recs in enumerate(it):
            if i % 50000 == 0:
                logger.info(
                    "Processed %d variants. Retained %d as they differ. Retained %d related to methylation.",
                    i,
                    len(this_res),
                    len(this_res_meth)
                )

            the_chrom = None
            the_ann = None
            for rec in recs:
                if rec is not None:
                    the_chrom = rec.CHROM
                    the_ann = '#'.join(rec.INFO['ANN'])
                    break

            if the_chrom in contigs:
                # compare across results
                gts = set()
                for rec in recs:
                    if rec is None:
                        this_res.append(dict(zip(readers.keys(), recs)))
                        # reset gts
                        gts = set()
                        break
                    else:
                        gts.add(rec.samples[0]['GT'])
                if len(gts) > 1:
                    this_res.append(dict(zip(readers.keys(), recs)))

                # regardless of whether they differ, if they are related to methylation keep track
                if any([t in the_ann for t in METH_GENES]):
                    this_res_meth.append(dict(zip(readers.keys(), recs)))

        out_fn = os.path.join(outdir, "%s_delta_variants.pkl.gz" % pid)
        with gzip.open(out_fn, 'wb') as f:
            pickle.dump(this_res, f)
        logger.info("Dumped %d results to gzipped pickle file %s", len(this_res), out_fn)

        out_fn = os.path.join(outdir, "%s_meth_variants.pkl.gz" % pid)
        with gzip.open(out_fn, 'wb') as f:
            pickle.dump(this_res_meth, f)
        logger.info("Dumped %d methylation-related results to gzipped pickle file %s", len(this_res_meth), out_fn)

