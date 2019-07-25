import math
import sys
import os
from pyscript import jobs, sge

try:
    PY_PATH = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'py')
except Exception:
    # probably running from an interpreter, so this won't work
    PY_PATH = '.'


class VCFRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 10


class VCFStatsBase(jobs.ArrayJob):
    title = 'vcf_stats'
    required_args = ['filename', 'variants_per_shard']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$FILE', '! -z'),
        ('$CONTIG', '! -z'),
        ('$STARTPOS', None),
        ('$ENDPOS', None),
    ]
    param_delim = ','

    core_cmd = """
        PATH="$PATH:%s"
        if [[ ! -z $STARTPOS ]]; then S="--start $STARTPOS"; else S=""; fi
        if [[ ! -z $ENDPOS ]]; then E="--end $ENDPOS"; else E=""; fi
        # get_gbm_specific_variants.py $FILE --contig $CONTIG $S $E --outdir {out_dir} {extra}
        get_methylation_related.py $FILE --contig $CONTIG $S $E --outdir {out_dir} {extra}
    """ % PY_PATH

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['filename'], self.args['variants_per_shard'])


class VCFStatsBash(jobs.BashArrayJobMixin, VCFStatsBase, jobs.VcfFileShardIterator):
    pass

class VCFStatsApocrita(VCFRequirements, VCFStatsBase, jobs.VcfFileShardIterator):
    pass


def run(run_type):
    import argparse
    import sys

    if run_type == 'bash':
        cls = VCFStatsBash
    elif run_type == 'apocrita':
        cls = VCFStatsApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    required.add_argument("filename", help="Path to VCF file")
    optional.add_argument("-n", "--variants_per_shard", help="Number of variants to handle per shard", type=int, default=int(1e7))
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-p", "--threads", help="Number of threads", default='1')

    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one here
        args.out_dir = os.path.join('./', 'vcf_stats')
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