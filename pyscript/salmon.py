import math
from pyscript import jobs, sge


class SalmonSgeRequirements(sge.ApocritaArrayJobMixin):
    @property
    def ram_per_core(self):
        eff_threads = int(self.args['threads']) + 1
        ram_per_core = 2.
        if (ram_per_core * eff_threads) < 8:
            ram_per_core = int(math.ceil(8. / float(eff_threads)))

        return "%dG" % ram_per_core

    @property
    def runtime_mins(self):
        return 60


class SalmonPEBase(jobs.ArrayJob):
    title = 'salmon_multilane_pe'
    required_args = ['read_dir', 'threads', 'index_dir', 'library_type']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
    ]
    param_delim = ':'
    core_cmd = "salmon quant -i {index_dir} -l {library_type} -p {threads} -1 $READS1 -2 $READS2 -o $SUBDIR"

    def set_default_extras(self):
        if '--seqBias' not in self.extra_args:
            self.extra_args.extend(
                ['--seqBias']
            )
        if '--gcBias' not in self.extra_args:
            self.extra_args.extend(
                ['--gcBias']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class SalmonSEBase(jobs.ArrayJob):
    title = 'salmon_multilane_se'
    required_args = ['read_dir', 'threads', 'index_dir', 'library_type']

    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$READ', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
    ]
    param_delim = ':'
    core_cmd = "salmon quant -i {index_dir} -l {library_type} -p {threads} -r $READ -o $SUBDIR"

    def set_default_extras(self):
        if '--seqBias' not in self.extra_args:
            self.extra_args.extend(
                ['--seqBias']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class SalmonMultilanePEApocrita(SalmonSgeRequirements, SalmonPEBase, jobs.PEFastqEncodeMultiLaneMixin):
    pass


class SalmonMultilanePEBash(jobs.BashArrayJobMixin, SalmonPEBase, jobs.PEFastqEncodeMultiLaneMixin):
    pass


class SalmonPEApocrita(SalmonSgeRequirements, SalmonPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class SalmonPEBash(jobs.BashArrayJobMixin, SalmonPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class SalmonSEApocrita(SalmonSgeRequirements, SalmonSEBase, jobs.SEFastqFileIteratorMixin):
    pass


class SalmonSEBash(jobs.BashArrayJobMixin, SalmonSEBase, jobs.SEFastqFileIteratorMixin):
    pass
