"""
Myrianni is designing gRNAs to target genes of interest that emerge from the computational pipeline. These will be
used together with a Crispr-CAS9/TET/DNMT system to (de)methylate key loci corresponding to the gene.

Aim here is to guide that process. Identify the genomic locations involved in any DMRs, and (ideally) the specific
probes that are most deregulated.
"""
import os
from matplotlib import pyplot as plt


from utils import output, genomics, log, setops
from plotting import genomics as genomic_plots
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from methylation import loader as methylation_loader, dmr, process


logger = log.get_console_logger()


if __name__ == '__main__':
    outdir = output.unique_output_dir()

    CELL_TYPE_MARKER_MAP = {
        'GBM': 'o',
        'iNSC': '^'
    }
    CELL_TYPE_COLOUR_MAP = {
        'GBM': 'k',
        'iNSC': 'w'
    }
    CELL_TYPE_ALPHA_MAP = {
        'GBM': 1,
        'iNSC': 0.4
    }
    CELL_TYPE_ZORDER_MAP = {
        'GBM': 19,
        'iNSC': 20
    }

    # load methylation and DMR data
    meth_obj = methylation_loader.load_by_patient(consts.PIDS, include_control=False)
    meth_obj.filter_by_sample_name(consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_INSC)
    meth_obj.meta.insert(0, 'patient_id', meth_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    mdat = process.m_from_beta(meth_obj.data)

    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')

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

    gois = ['PTGER4']
    tax_id = 9606
    genome_version = 'GRCh37'
    fa_fn = genomics.get_reference_genome_fasta(tax_id, version=genome_version)
    for g in gois:
        feat_res = genomics.get_features_from_gtf(g, tax_id=tax_id, version=genome_version)
        for i, (ftr, x) in enumerate(feat_res.items()):
            # get sequence
            logger.info("Gene %s, feature %d (%s). Strand: %s.", g, i + 1, ftr.id, ftr.strand)
            gene_seq = ftr.sequence(fa_fn)
            out_fn = os.path.join(outdir, "%s_full_sequence.txt" % g)
            with open(out_fn, 'wb') as f:
                f.write(gene_seq)
            logger.info("Wrote full genomic sequence (length %d) to %s.", len(gene_seq), out_fn)

            # get clusters and DMRs
            check_intersection = lambda x: \
                (ftr.start <= x.coord_range[0] <= ftr.stop) or \
                (ftr.start <= x.coord_range[1] <= ftr.stop) or \
                (x.coord_range[0] < ftr.start and x.coord_range[1] > ftr.stop)
            this_clusters = dict([
                (k, v) for k, v in dmr_res_s1.clusters.items() if v.chr == ftr.chrom and check_intersection(v)
            ])
            this_dmrs = {}
            for k, v in this_clusters.items():
                for pid in dmr_res_s1.iterkeys():
                    t = dmr_res_s1[pid].results[k]
                    if t.get('rej_h0', False):
                        this_dmrs.setdefault(k, {})[pid] = t

        # patients_involved = sorted(setops.reduce_union(*[t.keys() for t in this_dmrs.values()]))

        # nrows = len(patients_involved) + 1  # one ax for each patient plus one for the track
        nrows = len(consts.PIDS) + 1  # one ax for each patient plus one for the track
        gs = plt.GridSpec(nrows=nrows, ncols=1)
        fig = plt.figure(figsize=(8, 6))
        track_ax = fig.add_subplot(gs[-1])
        genomic_plots.plot_gene_transcripts(feat_res, ax=track_ax)
        # dm_axs = dict([(pid, fig.add_subplot(gs[i], sharex=track_ax)) for i, pid in enumerate(patients_involved)])
        dm_axs = dict([(pid, fig.add_subplot(gs[i], sharex=track_ax)) for i, pid in enumerate(consts.PIDS)])

        for cid, d1 in this_dmrs.items():
            pc = this_clusters[cid]
            for pid, d2 in d1.items():
                this_ax = dm_axs[pid]
                # bar showing DMR coverage and delta
                this_ax.barh(
                    0,
                    pc.coord_range[1] - pc.coord_range[0],
                    left=pc.coord_range[0],
                    height=d2['median_change'],
                    color=consts.PATIENT_COLOURS[pid],
                    align='edge',
                    zorder=15,
                    edgecolor='k',
                    linewidth=0.75
                )

        # scatter showing individual probe values
        ix = (anno.CHR == ftr.chrom) & (anno.MAPINFO >= ftr.start) & (anno.MAPINFO <= ftr.stop)
        this_probe_ids = anno.loc[anno.index[ix].intersection(mdat.index), 'MAPINFO'].sort_values().index
        for pid in consts.PIDS:
            this_ax = dm_axs[pid]
            for col, x in mdat.loc[this_probe_ids, meth_obj.meta.patient_id == pid].iteritems():
                the_typ = meth_obj.meta.type[col]
                this_ax.scatter(
                    anno.MAPINFO[this_probe_ids],
                    x.values,
                    c=CELL_TYPE_COLOUR_MAP[the_typ],
                    marker=CELL_TYPE_MARKER_MAP[the_typ],
                    zorder=CELL_TYPE_ZORDER_MAP[the_typ],
                    alpha=CELL_TYPE_ALPHA_MAP[the_typ]
                )




