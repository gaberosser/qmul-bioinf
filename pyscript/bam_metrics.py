import math
from pyscript import jobs, sge


class BamMetricsSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 30


class BamMetricsBase(jobs.ArrayJob):
    title = 'bam_metrics'
    required_args = ['bam_dir']

    ## FIXME: params file is more complex than this script expects
    # self.params.append([base, os.path.abspath(os.path.join(read_dir, t)), out_subdir])
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
        ('$BAM', '! -z'),
        ('$OUTFILE', '! -z'),
    ]

    ## FIXME: this doesn't remove the .sorted extension?

    core_cmd = """
    outfile="${{OUTFILE}}.bam_metrics"
    samtools view -b {extra} $BAM | samtools flagstat - > $outfile
    plus_counts=$(samtools view -F16 {extra} $BAM | cut -f3-4 | sort -u | wc -l)
    minus_counts=$(samtools view -f16 {extra} $BAM | cut -f3-4 | sort -u | wc -l)
    echo "Forward strand unique reads: $plus_counts" >> $outfile
    echo "Reverse strand unique reads: $minus_counts" >> $outfile
    total_counts=$(cat $outfile | grep 'in total (QC-passed' | cut -d ' ' -f1)
    NRF=`echo "($plus_counts + $minus_counts)/$total_counts" | bc -l`
    echo "NRF = $NRF" >> $outfile
    samtools stats {extra} $BAM >> $outfile
    """

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['bam_dir'])


class BamMetricsApocrita(BamMetricsSgeRequirements, BamMetricsBase, jobs.BamFileIteratorMixin):
    pass


class BamMetricsBash(jobs.BashArrayJobMixin, BamMetricsBase, jobs.BamFileIteratorMixin):
    pass


def multilane_run(run_type):
    import argparse
    import sys
    import os

    if run_type == 'bash':
        cls = BamMetricsBash
    elif run_type == 'apocrita':
        cls = BamMetricsApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--bam_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")

    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.bam_dir, 'bam_metrics')
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    if args.include is not None:
        args.include = args.include.split(',')
    if args.exclude is not None:
        args.exclude = args.exclude.split(',')

    obj = cls(extra_args=extra, **args.__dict__)
    obj.create_script()
    obj.submit()
