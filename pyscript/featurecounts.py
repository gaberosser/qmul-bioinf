from pyscript import jobs, sge


class FeatureCountsSgeRequirements(sge.ApocritaJobMixin):
    @property
    def ram_per_core(self):
        return "1G"

    @property
    def runtime_mins(self):
        return 60


class FeatureCountsBase(jobs.Job):
    title = 'featurecounts'
    required_args = ['read_dir', 'threads', 'gtf']

    param_delim = ':'
    core_cmd = "featureCounts -a {gtf} -T {threads} -o {out_dir}/quant {extra}"

    def core_command(self):
        cmd = super(FeatureCountsBase, self).core_command()
        return "{cmd} {bams}".format(cmd=cmd, bams=self.params)


    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])
        # force a single job (no array)
        self.params = self.params[1]


class FeatureCountsApocrita(FeatureCountsSgeRequirements, FeatureCountsBase, jobs.BamFileSingleRunMixin):
    pass


class FeatureCountsBash(jobs.BashJobMixin, FeatureCountsBase, jobs.BamFileSingleRunMixin):
    pass


def run(run_type):
    import argparse
    import sys
    import os

    if run_type == 'bash':
        cls = FeatureCountsBash
    elif run_type == 'apocrita':
        cls = FeatureCountsApocrita
    else:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-T", "--threads", help="Number of threads", default='1')
    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    required.add_argument("-a", "--gtf", help="GTF annotation file", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.read_dir, 'featurecounts')
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