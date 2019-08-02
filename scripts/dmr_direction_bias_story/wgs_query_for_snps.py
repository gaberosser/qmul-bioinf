import os
import gzip
import pickle
import vcf
import vcf.utils
import collections
import pandas as pd
from utils import setops, output, log, dictionary
from settings import DATA_DIR
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
        elif (t[key1].samples[0]['GT'] == '1/1') and (t[key2].samples[0]['GT'] == '0/1'):
            this_var_count[t[key1].CHROM]['%s hom %s het' % (key1, key2)] += 1
        elif (t[key2].samples[0]['GT'] == '1/1') and (t[key1].samples[0]['GT'] == '0/1'):
            this_var_count[t[key1].CHROM]['%s hom %s het' % (key2, key1)] += 1
        else:
            this_var_count[t[key1].CHROM]['other'] += 1
    return this_var_count


def run_one_sort(res_arr, key1, key2):
    res = collections.defaultdict(list)
    for t in res_arr:
        if t[key2] is None:
            res['%s only' % key1].append(t[key1])
        elif t[key1] is None:
            res['%s only' % key2].append(t[key2])
        elif (t[key1].samples[0]['GT'] == '1/1') and (t[key2].samples[0]['GT'] == '0/1'):
            res['%s hom %s het' % (key1, key2)].append(t)
        elif (t[key2].samples[0]['GT'] == '1/1') and (t[key1].samples[0]['GT'] == '0/1'):
            res['%s hom %s het' % (key2, key1)].append(t)
        else:
            res['other'].append(t)
    return res


if __name__ == '__main__':
    """
    Here we are going to iterate over a VCF file containing all variants called in our WGS dataset.
    We are looking specifically for any alterations in or in the vicinity of a pre-defined list of candidate genes
    known to be involved in DNA methylation.

    While we are at it, we may generate some overview plots, too.
    """
    pids = consts.PIDS
    contigs = set(['chr%d' % i for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])

    # V1: iterate over full VCF files and count (only - too big)
    if False:

        base_indir = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12')
        meta_fn = os.path.join(base_indir, 'sources.csv')

        meta = pd.read_csv(meta_fn, header=0, index_col=0)
        meta.loc[:, 'patient_id'] = ["%03d" % t for t in meta.patient_id]

        outdir = output.unique_output_dir()

        var_dat = []

        var_counts = collections.OrderedDict()
        meth_counts = collections.OrderedDict()

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

            var_counts[pid] = run_one_count(this_res)
            meth_counts[pid] = run_one_count(this_res_meth)

    # V2: iterate over pre-made short files and store data in memory
    base_indir = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12/meth_associated/')
    meta_fn = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12/', 'sources.csv')

    meta = pd.read_csv(meta_fn, header=0, index_col=0)
    meta.loc[:, 'patient_id'] = ["%03d" % t for t in meta.patient_id]

    outdir = output.unique_output_dir()

    var_dat = {}

    for pid in pids:
        logger.info("Patient %s", pid)
        this_meta = meta.loc[meta.patient_id == pid]
        in_fns = []
        readers = {}
        for t in this_meta.index:
            the_fn = os.path.join(base_indir, "%s.vcf" % meta.loc[t, 'sample'])
            in_fns.append(the_fn)
            readers[this_meta.loc[t, 'type']] = vcf.Reader(filename=the_fn)
        logger.info("Found %d conditions to compare: %s", len(readers), ', '.join(readers.keys()))

        it = vcf.utils.walk_together(*readers.values())

        this_res = []

        for i, recs in enumerate(it):
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

            var_dat[pid] = this_res

    dat_classified = dict([
        (pid, run_one_sort(var_dat[pid], 'GIC', 'iNSC')) for pid in pids
    ])

    # search through GIC only and GIC hom/iNSC het SNPs and 'other' and generate upset
    members = {}
    for pid, d in dat_classified.items():
        members[pid] = set()
        for typ in ['GIC only', 'GIC hom iNSC het', 'other']:
            for x in d[typ]:
                if isinstance(x, dict):
                    members[pid].add(str(x['GIC']))
                else:
                    members[pid].add(str(x))

    vs, vc = setops.venn_from_arrays(*[members[pid] for pid in pids])
    groups = {
        'Hypo': ['019', '030', '031', '017'],
        'Hyper': ['018', '050', '054', '061', '026', '052']
    }
    group_ind = setops.groups_to_ind(pids, groups)
    groups_inv = dictionary.complement_dictionary_of_iterables(groups, squeeze=True)

    venn_sets_by_group = setops.full_partial_unique_other_sets_from_groups(pids, groups)