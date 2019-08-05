import os
import gzip
import pickle
import itertools
import vcf.utils
import collections
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from utils import setops, output, log, dictionary
from settings import DATA_DIR
from scripts.hgic_final import consts


logger = log.get_console_logger()

METH_GENES = [
    "A0A096LPK6", "AICDA", "ALKBH1", "ALKBH2", "ALKBH3", "APEX1", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B",
    "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H", "ASZ1", "ATF7IP", "ATRX", "BAZ2A", "BEND3", "BRCA1",
    "CTCF", "CTCFL", "DDX4", "DMAP1", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DPPA3", "EHMT1", "EHMT2", "EZH2",
    "FAM129B", "FKBP6", "FOS", "FTO", "GATA3", "GATAD2A", "GNAS", "GRHL2", "GSK3A", "GSK3B", "H1FOO", "HELLS",
    "HEMK1",
    "KCNQ1OT1", "KDM1B", "KMT2A", "KMT2E", "MAEL", "MBD1", "MBD2", "MBD3", "MECP2", "METTL4", "MGMT", "MIS18A",
    "MORC1",
    "MOV10L1", "MPHOSPH8", "MTA2", "MTRR", "MYC", "N6AMT1", "OTUD4", "PARP1", "PICK1", "PIK3CA", "PIKC3A", "PIWIL2",
    "PIWIL4", "PLD6", "PPM1D", "PRDM14", "PRMT5", "PRMT7", "RLF", "SPI1", "STPG4", "TDG", "TDRD1", "TDRD12",
    "TDRD5",
    "TDRD9", "TDRKH", "TET1", "TET2", "TET3", "TRIM28", "UHRF1", "UHRF2", "USP7", "USP9X", "WT1", "ZFP57",
    "ZMPSTE24"
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

    hypo_count_full = vc[venn_sets_by_group['full']['Hypo'][0]]
    hyper_count_full = vc[venn_sets_by_group['full']['Hyper'][0]]

    hypo_counts_partial = [
        (setops.key_to_members(t, pids), vc[t]) for t in venn_sets_by_group['partial']['Hypo']
        ]
    hyper_counts_partial = [
        (setops.key_to_members(t, pids), vc[t]) for t in venn_sets_by_group['partial']['Hyper']
        ]

    # is this significant in any way?
    # focus on 3/4 of hypo OR 5/6 of hyper
    it = itertools.combinations(pids, 4)
    perm_hypo_full = []
    perm_hyper_full = []
    perm_hypo_partial = []
    perm_hyper_partial = []

    for fake_hypo_group in it:
        fake_groups = {
            'Hypo': fake_hypo_group,
            'Hyper': set(pids).difference(fake_hypo_group)
        }
        venn_sets_by_fake_group = setops.full_partial_unique_other_sets_from_groups(pids, fake_groups)

        perm_hypo_full.append(vc[venn_sets_by_fake_group['full']['Hypo'][0]])
        perm_hyper_full.append(vc[venn_sets_by_fake_group['full']['Hyper'][0]])

        fake_hypo_counts_partial = [
            (setops.key_to_members(t, pids), vc[t]) for t in venn_sets_by_fake_group['partial']['Hypo']
        ]
        fake_hyper_counts_partial = [
            (setops.key_to_members(t, pids), vc[t]) for t in venn_sets_by_fake_group['partial']['Hyper']
        ]

        perm_hypo_partial.append(
            sum([t[1] for t in fake_hypo_counts_partial if len(t[0]) > 2])
        )
        perm_hyper_partial.append(
            sum([t[1] for t in fake_hyper_counts_partial if len(t[0]) > 4])
        )

    fig, ax = plt.subplots()
    sns.kdeplot(perm_hypo_partial, ax=ax)
    ax.axvline(sum([t[1] for t in hypo_counts_partial if len(t[0]) > 2]), c='k', ls='--')
    ax.set_xlabel('Number of variants')
    ax.set_ylabel('Density')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "permute_partial_counts_meth_assoc_hypo.png"), dpi=200)

    fig, ax = plt.subplots()
    sns.kdeplot(perm_hyper_partial, ax=ax)
    ax.axvline(sum([t[1] for t in hyper_counts_partial if len(t[0]) > 4]), c='k', ls='--')
    ax.set_xlabel('Number of variants')
    ax.set_ylabel('Density')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "permute_partial_counts_meth_assoc_hyper.png"), dpi=200)

    # track these partial matches down
    aa_hypo = [(setops.key_to_members(t, pids)[0], vs[t]) for t in venn_sets_by_group['partial']['Hypo'] if len(setops.key_to_members(t, pids)) > 2]
    aa_hyper = [(setops.key_to_members(t, pids)[0], vs[t]) for t in venn_sets_by_group['partial']['Hyper'] if len(setops.key_to_members(t, pids)) > 4]
    all_hypo = setops.reduce_union(*[t[1] for t in aa_hypo])
    all_hyper = setops.reduce_union(*[t[1] for t in aa_hyper])

    partial_hypo_recs = []
    partial_hyper_recs = []

    for pid, arr in aa_hypo:
        for x in arr:
            the_search_list = dat_classified[pid]['GIC only']
            the_search_list.extend([t['GIC'] for t in dat_classified[pid]['GIC hom iNSC het']])
            the_recs = [t for t in the_search_list if str(t) == x]
            if len(the_recs) > 0:
                partial_hypo_recs.append(the_recs[0])
            else:
                print "No search result for rec %s" % x

    for pid, arr in aa_hyper:
        for x in arr:
            the_search_list = dat_classified[pid]['GIC only']
            the_search_list.extend([t['GIC'] for t in dat_classified[pid]['GIC hom iNSC het']])
            the_recs = [t for t in the_search_list if str(t) == x]
            if len(the_recs) > 0:
                partial_hyper_recs.append(the_recs[0])
            else:
                print "No search result for rec %s" % x

    # V3: Iterate over ALL VCFs in one and look for fully group-specific variants
    group_hypo = set(groups['Hypo'])
    group_hyper = set(groups['Hyper'])

    base_indir = os.path.join(DATA_DIR, 'wgs', 'x17067/2017-12-12')
    gs_var = collections.defaultdict(list)
    gs_var_hyper = []
    readers = {}

    for pid in pids:
        # logger.info("Patient %s", pid)
        this_meta = meta.loc[meta.patient_id == pid]
        for t in this_meta.index:
            the_fn = os.path.join(base_indir, t, "%s.vcf.gz" % meta.loc[t, 'sample'])
            readers[(pid, this_meta.loc[t, 'type'])] = vcf.Reader(filename=the_fn, compressed=True)

    it = vcf.utils.walk_together(*readers.values())
    count = 0

    for recs in it:
        count += 1

        if count % 50000 == 0:
            logger.info(
                "Processed %d variants.",
                count,
            )

        rec_dict = dict(zip(readers.keys(), recs))
        the_chrom = None
        the_ann = None
        for rec in recs:
            if rec is not None:
                the_chrom = rec.CHROM
                the_ann = '#'.join(rec.INFO['ANN'])
                break

        if the_chrom not in contigs:
            continue

        this_pids = []

        for pid in pids:
            a = rec_dict[(pid, 'GIC')]
            b = rec_dict[(pid, 'iNSC')]
            if a is None:
                continue
            elif b is None:
                this_pids.append(pid)
            elif (a.samples[0]['GT'] == '1/1') and (b.samples[0]['GT'] == '0/1'):
                this_pids.append(pid)

        if set(this_pids) == group_hypo:
            gs_var['hypo'].append(rec_dict)
            logger.info("Found %d hypo-related variants in %d records", len(gs_var['hypo']), count)
        elif set(this_pids) == group_hyper:
            gs_var['hyper'].append(rec_dict)
            logger.info("Found %d hyper-related variants in %d records", len(gs_var['hyper']), count)

    # TEMP
    with open(os.path.join('/home/gabriel/Dropbox/pigeonhole', "all_variants_group_specific.pkl"), 'wb') as f:
        pickle.dump(gs_var, f)

    with open(os.path.join(outdir, "all_variants_group_specific.pkl"), 'wb') as f:
        pickle.dump(gs_var, f)