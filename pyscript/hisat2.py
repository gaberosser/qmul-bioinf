import math
from pyscript import jobs, sge


class Hisat2SgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "%dG" % 4

    @property
    def runtime_mins(self):
        return 240


class Hisat2SEBase(jobs.ArrayJob):
    title = 'hisat2_se'
    required_args = ['read_dir', 'threads', 'index_dir']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$NAME', None),
        ('$READ', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
    ]
    param_delim = ':'
    core_cmd = "hisat2 -x {index_dir} -p {threads} -U $READ -S $SUBDIR"

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class Hisat2SEApocrita(Hisat2SgeRequirements, Hisat2SEBase, jobs.SEFastqFileIteratorMixin):
    pass


class Hisat2SEBash(jobs.BashArrayJobMixin, Hisat2SEBase, jobs.SEFastqFileIteratorMixin):
    pass
