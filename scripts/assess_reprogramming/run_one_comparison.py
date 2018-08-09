import sys
import os
import pickle
import multiprocessing as mp
# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')

import pandas as pd
import argparse

DMR_PARAMS = {
    'd_max': 400,
    'n_min': 6,
    'delta_m_min': 1.4,
    'alpha': 0.01,
    'dmr_test_method': 'mwu',  # 'mwu', 'mwu_permute'
    'test_kwargs': {},
    'n_jobs': mp.cpu_count(),
}


def run_dmr(dmr_clusters, dat, samples1, samples2, **dmr_params):
    """
    Run DMR
    :param dmr_params:
    :return:
    """
    the_samples = [samples1, samples2]
    dmr_clusters.test_clusters(
        dat,
        samples=the_samples,
        n_jobs=dmr_params['n_jobs'],
        min_median_change=dmr_params['delta_m_min'],
        method=dmr_params['dmr_test_method'],
        alpha=dmr_params['alpha'],
        **dmr_params['test_kwargs']
    )
    return dmr_clusters


def main():
    parser = argparse.ArgumentParser(description="Run a single DMR comparison")
    parser.add_argument("--clusters", help="path to clusters pkl file", required=True)
    parser.add_argument("--name", help="name for this comparison", required=True)
    parser.add_argument("-t","--threads", help="number of threads", default=1, type=int)
    parser.add_argument("-O","--outdir", help="output directory", default=".")
    parser.add_argument("--data_dir", help="directory containing methylation data files", required=True)
    parser.add_argument("--samples_1", help="comma-separated list of sample names for group 1", required=True)
    parser.add_argument("--samples_2", help="comma-separated list of sample names for group 2", required=True)
    args, extra = parser.parse_known_args()

    if not os.path.isdir(args.outdir):
        sys.stderr.write("Requested output dir %s does not exist, creating.\n" % args.outdir)
        os.makedirs(args.outdir)

    if not os.path.isfile(args.clusters):
        sys.stderr.write("Invalid clusters file: %s\n" % args.clusters)
        sys.exit(1)

    if not os.path.isdir(args.data_dir):
        sys.stderr.write("Invalid data directory: %s\n" % args.data_dir)
        sys.exit(1)

    with open(args.clusters, 'rb') as f:
        dmr_clusters = pickle.load(f)
    sys.stderr.write("Loaded %d DMR clusters to be tested.\n" % len(dmr_clusters.clusters))

    samples_1 = args.samples_1.split(',')
    samples_2 = args.samples_2.split(',')

    sys.stderr.write("Requested a %d vs %d DMR comparison.\n" % (len(samples_1), len(samples_2)))
    sys.stderr.write("Group 1: %s\n" % args.samples_1)
    sys.stderr.write("Group 2: %s\n" % args.samples_2)

    # load sample data
    dat = None
    for s in samples_1 + samples_2:
        fn = os.path.join(args.data_dir, "%s.csv.gz" % s)
        if not os.path.isfile(fn):
            sys.stderr.write("ERROR: No file found at %s.\n" % fn)
            sys.exit(1)
        this = pd.read_csv(fn, header=None, index_col=0)
        this.index.name = s
        if dat is None:
            dat = this
            dat.columns = [s]
        else:
            dat.insert(dat.shape[1], s, this)

    sys.stderr.write("Successfully loaded %d samples with %d probes.\n" % (dat.shape[1], dat.shape[0]))
    res = run_dmr(dmr_clusters, dat, samples_1, samples_2, **DMR_PARAMS)
    sys.stderr.write("Successfully run DMR, resulting in %d significant clusters.\n" % len(res.results_significant))
    tbl = res.to_table(skip_geneless=False)

    fn_out = os.path.join(args.outdir, "%s.csv" % args.name)
    tbl.to_csv(fn_out)
    sys.stderr.write("Wrote results table to output file %s.\n" % fn_out)
    sys.exit(0)


if __name__ == "__main__":
    main()
