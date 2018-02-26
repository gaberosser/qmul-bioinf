import math
from pyscript import jobs, sge

class BismarkSgeRequirements(sge.ApocritaArrayJobMixin):
    pass
    ## TODO
    # min_ram_gb = 32.
    # estimated_runtime_per_core = 720.
    # @property
    # def ram_per_core(self):
    #     nthread = int(self.args['threads'])
    #     # aim for minimum 32Gb
    #     gb_per_core = int(math.ceil(self.min_ram_gb / float(nthread)))
    #     print "Num thread: %d" % nthread
    #     print "Gb per core: %d" % gb_per_core
    #     return "%dG" % gb_per_core
    #
    # @property
    # def runtime_mins(self):
    #     nthread = int(self.args['threads'])
    #     # roughly 60 mins with 12 cores
    #     return self.estimated_runtime_per_core / float(nthread)


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


class BismarkMultilanePEBash(jobs.BashArrayJobMixin, BismarkPEBase, jobs.PEFastqBartsMultiLaneMixin):
    pass