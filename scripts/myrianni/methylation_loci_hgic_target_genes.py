"""
Myrianni is designing gRNAs to target genes of interest that emerge from the computational pipeline. These will be
used together with a Crispr-CAS9/TET/DNMT system to (de)methylate key loci corresponding to the gene.

Aim here is to guide that process. Identify the genomic locations involved in any DMRs, and (ideally) the specific
probes that are most deregulated.
"""
import os
from matplotlib import pyplot as plt
import pysam
import pandas as pd
import pickle

from utils import output, genomics, log, setops
from plotting import genomics as genomic_plots
from scripts.hgic_final import consts, two_strategies_grouped_dispersion as tsgd
from methylation import loader as methylation_loader, dmr, process
from rnaseq import loader as rnaseq_loader


logger = log.get_console_logger()


if __name__ == '__main__':
    outdir = output.unique_output_dir()
    edge_buffer_bp = 1500

    CELL_TYPE_MARKER_MAP = {
        'GBM': 'o',
        'iNSC': 'o'
    }
    CELL_TYPE_COLOUR_MAP = {
        'GBM': 'r',
        'iNSC': 'b'
    }
    CELL_TYPE_ALPHA_MAP = {
        'GBM': 0.7,
        'iNSC': 0.6
    }
    CELL_TYPE_ZORDER_MAP = {
        'GBM': 19,
        'iNSC': 20
    }
    CELL_TYPE_SIZE_MAP = {
        'GBM': 25,
        'iNSC': 20
    }

    # load methylation and DMR data
    meth_obj = methylation_loader.load_by_patient(consts.PIDS, include_control=False)
    meth_obj.filter_by_sample_name(consts.S1_METHYL_SAMPLES_GIC + consts.S1_METHYL_SAMPLES_INSC)
    meth_obj.meta.insert(0, 'patient_id', meth_obj.meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))

    mdat = process.m_from_beta(meth_obj.data)

    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS
    de_params = consts.DE_PARAMS

    DMR_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'dmr')
    DE_LOAD_DIR = os.path.join(output.OUTPUT_DIR, 'de')

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

    gois = ['PTGER4', 'NTRK2']
    rois = {
        'PTGER4': [[40679018, 40682333]],  # dmr_res_s1.clusters[5877].coord_range
        'NTRK2': [[87283670, 87287300]]
    }
    tax_id = 9606
    genome_version = 'GRCh37'
    fa_fn = genomics.get_reference_genome_fasta(tax_id, version=genome_version, ext=['.fa', '.fa.bgz'])
    for g in gois:
        this_clusters = dict([
            (k, v) for k, v in dmr_res_s1.clusters.items() if len([x for x in v.genes if x[0] == g])
        ])

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
            this_clusters.update(
                dict([
                    (k, v) for k, v in dmr_res_s1.clusters.items() if v.chr == ftr.chrom and check_intersection(v)
                ])
            )
        this_dmrs = {}
        for k, v in this_clusters.items():
            for pid in dmr_res_s1.iterkeys():
                t = dmr_res_s1[pid].results[k]
                if t.get('rej_h0', False):
                    this_dmrs.setdefault(k, {})[pid] = t

        nrows = len(consts.PIDS) + 1  # one ax for each patient plus one for the track
        # first column shows methylation data, second RNA-Seq DE
        gs = plt.GridSpec(nrows=nrows, ncols=2, width_ratios=(10, 1))
        fig = plt.figure(figsize=(8, 6))
        track_ax = fig.add_subplot(gs[-1, 0])
        sharex = None
        de_share = None
        genomic_plots.plot_gene_transcripts(feat_res, ax=track_ax)
        dm_axs = {}
        de_axs = {}
        for i, pid in enumerate(consts.PIDS):
            dm_axs[pid] = fig.add_subplot(gs[i, 0], sharex=sharex)
            sharex = dm_axs[pid]
            de_axs[pid] = fig.add_subplot(gs[i, 1], sharex=de_share, sharey=de_share)
            de_share = de_axs[pid]

        # keep tabs on all probes involved in DMRs (they may be outside of the region)
        probes_in_dmrs = set()

        for cid, d1 in this_dmrs.items():
            pc = this_clusters[cid]
            probes_in_dmrs.update(pc.pids)
            for pid, d2 in d1.items():
                this_ax = dm_axs[pid]
                # bar showing DMR coverage and delta
                this_ax.barh(
                    0,
                    pc.coord_range[1] - pc.coord_range[0],
                    left=pc.coord_range[0],
                    height=d2['median_change'],
                    color=consts.METHYLATION_DIRECTION_COLOURS['hyper' if d2['median_change'] > 0 else 'hypo'],
                    align='edge',
                    zorder=15,
                    edgecolor='k',
                    linewidth=0.75
                )

        # scatter showing individual probe values
        ix = (
            (anno.CHR == ftr.chrom) & (anno.MAPINFO >= ftr.start) & (anno.MAPINFO <= ftr.stop)
        ) | (anno.index.isin(probes_in_dmrs))
        this_probe_ids = anno.loc[anno.index[ix].intersection(mdat.index), 'MAPINFO'].sort_values().index
        ymin = 0
        ymax = 0
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
                    alpha=CELL_TYPE_ALPHA_MAP[the_typ],
                    s=CELL_TYPE_SIZE_MAP[the_typ],
                    edgecolor='k',
                    linewidth=0.5
                )
                ymin = min(x.values.min(), ymin)
                ymax = max(x.values.max(), ymax)
                this_ax.set_ylabel(pid)

            this_de_ax = de_axs[pid]
            # bar chart to represent DE
            ix = de_res_s1[pid]['Gene Symbol'] == g
            if ix.sum() == 1:
                this_de = de_res_s1[pid].loc[ix].squeeze()
                this_de_ax.barh(
                    [0],
                    [this_de.logFC],
                    edgecolor='k' if this_de.FDR < de_params['fdr'] else 'none',
                    linewidth=1.,
                    color=consts.METHYLATION_DIRECTION_COLOURS['hyper'] if this_de.logFC > 0 else consts.METHYLATION_DIRECTION_COLOURS['hypo']
                )
            elif ix.sum() > 1:
                logger.warn("Gene %s. Found %d matching DE results in patient %s??", g, int(ix.sum()), pid)
            else:
                logger.warn("Gene %s. No matching DE result for patient %s.", g, pid)


        dm_axs[consts.PIDS[-1]].xaxis.set_ticklabels([])
        plt.setp([t.yaxis for t in de_axs.values()], visible=False)
        de_axs[consts.PIDS[-1]].set_ylim([-1, 1])
        plt.setp(de_axs[consts.PIDS[-1]].xaxis.get_ticklabels(), rotation=90)
        de_axs[consts.PIDS[-1]].set_xlabel('DE logFC')

        [ax.set_ylim([ymin * 1.05, ymax * 1.05]) for ax in dm_axs.values()]

        gs.update(top=0.98, left=0.07, right=0.97, bottom=0.09, hspace=0.15, wspace=0.05)
        fig.savefig(os.path.join(outdir, "%s_mvalues.png" % g), dpi=300)

        if g in rois:
            for i, (x0, x1) in enumerate(rois[g]):
                track_ax.set_xlim([x0 - 200, x1 + 200])
                dm_axs.values()[0].set_xlim([x0 - 200, x1 + 200])
                fig.savefig(os.path.join(outdir, "%s_mvalues_roi_%d.png" % (g, i + 1)), dpi=300)



        # by request, create a file containing coord, seq and indicating the probe locations
        the_ftr = [k for k in feat_res.keys() if k.featuretype == 'gene']
        if len(the_ftr) != 1:
            raise NotImplementedError("We require exactly one gene feature for this next step")
        the_ftr = the_ftr[0]
        # expand and take the sequence
        ff = pysam.Fastafile(fa_fn)
        the_seq = ff.fetch('5', the_ftr.start - edge_buffer_bp - 1, the_ftr.stop + edge_buffer_bp)
        # the_seq = the_ftr.sequence(fa_fn)

        the_coords = range(the_ftr.start - edge_buffer_bp, the_ftr.stop + edge_buffer_bp + 1)
        if the_ftr.strand == '-':
            logger.warn("FIXME: wheck the implementation of reverse-stranded genes")
            the_coords = the_coords[::-1]
        ix = (
            (anno.CHR == the_ftr.chrom) & (anno.MAPINFO >= the_ftr.start) & (anno.MAPINFO <= the_ftr.stop)
        ) | (anno.index.isin(probes_in_dmrs))
        this_probe_ids = anno.loc[anno.index[ix].intersection(mdat.index), 'MAPINFO'].sort_values().index
        this_probe_coords = anno.loc[this_probe_ids, 'MAPINFO']
        this_probe_indicator = [0] * len(the_coords)
        for i, cc in enumerate(the_coords):
            if cc in this_probe_coords.values:
                this_probe_indicator[i] = 1

        # create a pd DataFrame for export
        for_export = pd.DataFrame(index=the_coords, dtype=str)
        for_export.insert(0, 'sequence', [t for t in the_seq])
        for_export.insert(1, 'probe', this_probe_indicator)
        i = 2
        for pid in consts.PIDS:
            this_mdat = mdat.loc[this_probe_ids, meth_obj.meta.patient_id == pid]
            for col, vals in this_mdat.iteritems():
                t = pd.Series(index=the_coords)
                t.loc[this_probe_coords.values] = vals.values
                for_export.insert(i, col, t)
                i += 1

        for_export.to_excel(os.path.join(outdir, "%s_sequence_coords_probes_values.xlsx" % g))
