import os
from matplotlib import pyplot as plt
import pickle
import collections

from utils import output, log, setops
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from methylation import loader as methylation_loader, dmr, process
from rnaseq import loader as rnaseq_loader
from settings import INTERMEDIATE_DIR

from plotting import genomics

logger = log.get_console_logger()


if __name__ == '__main__':
    outdir = output.unique_output_dir()

    # load methylation and DMR data
    meth_obj = methylation_loader.load_by_patient(consts.PIDS, include_control=False)
    meth_obj.filter_by_sample_name(consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_INSC)
    meth_obj.meta.insert(0, 'patient_id', meth_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    mdat = process.m_from_beta(meth_obj.data)

    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    de_params = consts.DE_PARAMS

    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'de')

    anno = methylation_loader.load_illumina_methylationepic_annotation()

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    # load DMR results
    the_hash = tsgd.dmr_results_hash(meth_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise Exception("Unable to locate pre-existing results.")

    # load DE results
    rnaseq_obj = rnaseq_loader.load_by_patient(consts.PIDS, include_control=False)
    rnaseq_obj.filter_by_sample_name(consts.S1_RNASEQ_SAMPLES)

    the_hash = tsgd.de_results_hash(rnaseq_obj.meta.index.tolist(), de_params)
    filename = 'de_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DE_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Reading S1 DE results from %s", fn)
        with open(fn, 'rb') as f:
            de_res_s1 = pickle.load(f)
    else:
        raise Exception("Unable to locate pre-existing DE results.")

    the_hash = tsgd.dmr_results_hash(meth_obj.meta.index.tolist(), dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
    else:
        raise Exception("Unable to locate pre-existing DMR results.")

    # prepare the DMR comparison groups dictionary required to specify which samples are being compared
    dmr_comparison_groups = collections.OrderedDict([(pid, {}) for pid in consts.PIDS])
    gg = mdat.columns.groupby(zip(meth_obj.meta.patient_id, meth_obj.meta.type))
    for (pid, typ), samples in gg.items():
        dmr_comparison_groups[pid][typ] = samples

    obj = genomics.MethylationExpressionLocusPlotter()
    obj.set_mvalues(mdat)
    obj.set_dmr_res(dmr_res_s1, dmr_comparison_groups)
    obj.set_de_res(de_res_s1)

    obj.plot_gene('PTGER4')