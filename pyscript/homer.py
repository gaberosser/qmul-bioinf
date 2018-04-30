import math
import csv
import os
import re
from pyscript import jobs, sge


class HomerAnnotatePeaksSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "4G"

    @property
    def runtime_mins(self):
        return 60


class HomerAnnotatePeaksBase(jobs.ArrayJob, jobs.Macs2PeaksIteratorMixin):
    title = 'homer_annotatepeaks'
    required_args = ['in_dir', 'fa_file']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
        ('$BED', '! -z'),
        ('$OUTSTEM', '! -z'),
    ]
    param_delim = ':'
    # TODO
    core_cmd = """
    annotatePeaks.pl $BED {fa_file} {extra} > ${{OUTSTEM}}.annotatePeaks
    """

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['in_dir'])


class HomerAnnotatePeaksApocrita(HomerAnnotatePeaksSgeRequirements, HomerAnnotatePeaksBase):
    core_cmd = """
    cp $BED $TMPDIR
    BED="$TMPDIR/$(basename $BED)"
    """ + HomerAnnotatePeaksBase.core_cmd


class HomerAnnotatePeaksBash(jobs.BashArrayJobMixin, HomerAnnotatePeaksBase):
    pass


def homer_annotate_peaks_run(run_type):
    import argparse
    import sys

    if run_type == 'bash':
        cls = HomerAnnotatePeaksBash
    elif run_type == 'apocrita':
        cls = HomerAnnotatePeaksApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("--in_dir", help="Path to peak files", default='./')
    parser.add_argument("--fa_file", help="Path to fasta reference", required=True)
    parser.add_argument("--include", help="Comma-separated list of filestems to include", required=False)
    parser.add_argument("--exclude", help="Comma-separated list of filestems to exclude", required=False)


    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.in_dir, 'homer_annotatepeaks')
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
