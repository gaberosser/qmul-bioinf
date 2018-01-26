import math
from pyscript import jobs, sge


class CufflinksSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        return "4G"

    @property
    def runtime_mins(self):
        return 360


class CufflinksBase(jobs.ArrayJob, jobs.BamFileIteratorMixin):
    title = 'cufflinks'
    required_args = ['bam_dir', 'threads', 'library_type', 'gtf']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$NAME', None),
        ('$BAM', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
    ]
    param_delim = ':'

    core_cmd = "cufflinks -G {gtf} -p {threads} -o $SUBDIR --library-type {library_type} $BAM"

    def set_default_extras(self):
        if '--no-update-check' not in self.extra_args:
            self.extra_args.extend(
                ['--no_update_check']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['bam_dir'])


class CufflinksApocrita(CufflinksSgeRequirements, CufflinksBase):
    pass


class CufflinksBash(jobs.BashArrayJobMixin, CufflinksBase):
    pass