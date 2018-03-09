import math
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
        return "{cmd} {bams}".format(cmd=cmd, bams=self.params[1])


    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class FeatureCountsApocrita(FeatureCountsSgeRequirements, FeatureCountsBase, jobs.BamFileSingleRunMixin):
    pass


class FeatureCountsBash(jobs.BashJobMixin, FeatureCountsBase, jobs.BamFileSingleRunMixin):
    pass
