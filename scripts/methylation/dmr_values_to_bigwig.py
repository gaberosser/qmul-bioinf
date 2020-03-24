import pyBigWig
import os
import collections
from methylation import loader, dmr, process
from scripts.hgic_final import two_strategies_grouped_dispersion as tsgd, consts
from utils import output, setops, genomics, log
from settings import LOCAL_DATA_DIR, INTERMEDIATE_DIR


logger = log.get_console_logger()


def write_bigwig(
        dmr_results,
        clusters,
        chrom_length,
        bw_fn,
        attr_key='median_change',
        chr_prefix=True
):
    bw = pyBigWig.open(bw_fn, 'w')
    bw.addHeader(chrom_length.items())

    # add entries, one chromosome at a time
    this_res = collections.defaultdict(list)
    for cid, attr in dmr_results.items():
        the_cluster = clusters[cid]
        the_chr = the_cluster.chr
        if chr_prefix:
            the_chr = "chr" + the_chr
        this_res[the_chr].append(the_cluster.coord_range + (attr[attr_key],))

    # order by genomic coordinate
    for chrom in chrom_length:
        if chrom not in this_res:
            continue
        x = sorted(this_res[chrom], key=lambda x: x[0])
        this_res[chrom] = x
        x_exp = zip(*x)

        # add to bigwig
        bw.addEntries([chrom] * len(x), list(x_exp[0]), ends=list(x_exp[1]), values=list(x_exp[2]))

    bw.close()


if __name__ == '__main__':

    """
    This script will create bigwig files for visualisation in a genome browser.
    These will contain the locations of DMRs (in GRCh37 coords) and the approximate DM median value
    Can also consider adding individual probes??
    """
    # if True, add 'chr' to all chromosome names - useful for cpompatibility with UCSC annotation track in IGV
    chr_prefix = True

    outdir = output.unique_output_dir()
    pids = consts.PIDS
    norm_method_s1 = 'swan'
    dmr_params = consts.DMR_PARAMS

    # set this to True if output bed files are required (this is quite slow due to the large number of combinations)
    write_bed_files = False

    outdir = output.unique_output_dir()
    DMR_LOAD_DIR = os.path.join(INTERMEDIATE_DIR, 'dmr')

    if False:
        me_obj, anno = tsgd.load_methylation(pids, norm_method=norm_method_s1, patient_samples=consts.S1_METHYL_SAMPLES)
        me_data = me_obj.data
        me_meta = me_obj.meta
        me_meta.insert(0, 'patient_id', me_meta.index.str.replace(r'(GBM|DURA)(?P<pid>[0-9]{3}).*', '\g<pid>'))
    else:
        anno = loader.load_illumina_methylationepic_annotation(split_genes=True)

    # We load pre-computed results if a file with the correct filename is found
    # Otherwise raise an error: this should be computed elsewhere

    # use a hash on the PIDs and parameters to ensure we're looking for the right results
    dmr_hash_dict = dict(dmr_params)
    dmr_hash_dict['norm_method'] = norm_method_s1

    the_hash = tsgd.dmr_results_hash(consts.S1_METHYL_SAMPLES, dmr_hash_dict)
    filename = 'dmr_results_paired_comparison.%d.pkl' % the_hash
    fn = os.path.join(DMR_LOAD_DIR, filename)

    if os.path.isfile(fn):
        logger.info("Loading pre-computed DMR results from %s", fn)
        dmr_res_s1 = dmr.DmrResultCollection.from_pickle(fn, anno=anno)
        clusters = dmr_res_s1.clusters
    else:
        raise IOError("Unable to locate pre-existing results.")

    # chromosome lengths
    chroms = [str(t) for t in range(1, 23)]
    fa_fn = os.path.join(
        LOCAL_DATA_DIR,
        'reference_genomes',
        'human/ensembl/GRCh38.release87/fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    )
    cl = genomics.feature_lengths_from_fasta(fa_fn, features=chroms)

    chrom_length = collections.OrderedDict()
    for k in chroms:
        if chr_prefix:
            new_k = "chr%s" % k
        else:
            new_k = k
        chrom_length[new_k] = cl[k]

    for pid in pids:
        fn = os.path.join(outdir, "%s_dmrs.bw" % pid)
        write_bigwig(
            dmr_res_s1[pid].results_significant,
            clusters,
            chrom_length,
            fn,
            chr_prefix=chr_prefix
        )

    # repeat for patient-specific DMRs
    patient_specific_cids = dict(zip(
        pids,
        setops.specific_features(*[dmr_res_s1[pid].results_significant for pid in pids])
    ))

    for pid in pids:
        fn = os.path.join(outdir, "%s_specific_dmrs.bw" % pid)
        this_res = dict([
            (cid, dmr_res_s1[pid].results[cid]) for cid in patient_specific_cids[pid]
        ])
        write_bigwig(
            this_res,
            clusters,
            chrom_length,
            fn,
            chr_prefix=chr_prefix
        )



