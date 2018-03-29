import math
import re
import os
from pyscript import jobs, sge


class SRASgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 180


class SRAGetterBase(jobs.ArrayJob):
    title = 'sra_getter'
    required_args = ['project_id', 'threads']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', None),
    ]
    param_delim = ':'
    core_cmd = "fastq-dump {extra} $ID -O {out_dir}"

    def set_default_extras(self):
        if '--gzip' not in self.extra_args:
            self.extra_args.extend(
                ['--gzip']
            )
        if '--split-files' not in self.extra_args:
            self.extra_args.extend(
                ['--split-files']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['project_id'])


class SRAGetterApocrita(SRASgeRequirements, SRAGetterBase, jobs.SraRunIteratorMixin):
    pass


class SRAGetterBash(jobs.BashArrayJobMixin, SRAGetterBase, jobs.SraRunIteratorMixin):
    pass


def run(run_type):
    import argparse
    import sys
    import os

    if run_type == 'bash':
        cls = SRAGetterBash
    elif run_type == 'apocrita':
        cls = SRAGetterApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    required.add_argument("--project_id", help="Project ID, PRJNA or SRP", required=True)
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-p", "--threads", help="Number of threads", default='1')

    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join('./', 'fastq')
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