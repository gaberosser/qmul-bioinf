import math
import sys
import os
from pyscript import jobs, sge

try:
    PY_PATH = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'py')
except Exception:
    # probably running from an interpreter, so this won't work
    PY_PATH = '.'


class MergeBamRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "2G"

    @property
    def runtime_mins(self):
        return 20


class MergeBamsBase(jobs.ArrayJob):
    title = 'merge_bams'
    required_args = ['read_dir', 'threads']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$BAMS', '! -z'),
        ('$OUTDIR', '! -z'),
    ]
    param_delim = ':'

    core_cmd = """
    outfn="$OUTDIR/${{ID}}.bam"
        samtools merge {extra} $outfn $BAMS -@ {threads}
    """

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class MergeBamsMultilaneBartsBash(jobs.BashArrayJobMixin, MergeBamsBase, jobs.BamBartsMultilaneRecursiveIteratorMixin):
    pass

class MergeBamsMultilaneBartsApocrita(MergeBamRequirements, MergeBamsBase, jobs.BamBartsMultilaneRecursiveIteratorMixin):
    pass


def run(run_type):
    import argparse
    import sys

    if run_type == 'bash_multilane_barts':
        cls = MergeBamsMultilaneBartsBash
    elif run_type == 'apocrita_multilane_barts':
        cls = MergeBamsMultilaneBartsApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    required.add_argument("read_dir", help="Input directory")
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("--threads", help="Number of threads to use", type=int, default=1)

    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one here
        args.out_dir = '.'
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    if args.include is not None:
        args.include = args.include.split(',')
    if args.exclude is not None:
        args.exclude = args.exclude.split(',')

    obj = cls(extra_args=extra, **args.__dict__)
    obj.create_script()
    obj.submit()