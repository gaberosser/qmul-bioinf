import math
from pyscript import jobs, sge

class BismarkSgeRequirements(sge.ApocritaArrayJobMixin):

    def shebang(self):
        """
        Generate the shebang line(s) for the top of the script
        :return:
        """
        # Default behaviour: run with bash
        return sge.apocrita_submission_header(
            work_dir=self.out_dir,
            threads=2 * int(self.args['threads']),
            ram_per_core=self.ram_per_core,
            runtime_mins=self.runtime_mins,
            arr_size=len(self.params)
        )

    @property
    def ram_per_core(self):
        return "2G"

    @property
    def runtime_mins(self):
        nthread = int(self.args['threads'])
        # roughly 30 mins with 4 threads
        estimated_runtime_one_core = 120.
        # allow for some kind of plateau
        return estimated_runtime_one_core / float(nthread) ** 0.8


class BismarkPEBase(jobs.ArrayJob):
    
    title = 'bismark'
    required_args = ['read_dir', 'threads', 'index_dir']
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
        ('$GZ', None)
    ]
    param_delim = ':'

    core_cmd = "bismark {extra} {index_dir} " + \
    "-o $SUBDIR -p {threads} -B $ID -1 $READS1 -2 $READS2"

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class BismarkPEBash(jobs.BashArrayJobMixin, BismarkPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class BismarkMultilanePEBash(jobs.BashArrayJobMixin, BismarkPEBase, jobs.PEFastqBartsMultiLaneMixin):
    ## FIXME: not sure this is working as it should
    pass


class BismarkPEApocrita(BismarkSgeRequirements, BismarkPEBase, jobs.PEFastqFileIteratorMixin):
    pass