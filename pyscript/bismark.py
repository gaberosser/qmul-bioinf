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

    ## FIXME: how do we pass extra args into two different calls?
    core_cmd = """
    bismark {extra} {index_dir} -o $SUBDIR -p {threads} -B $ID -1 $READS1 -2 $READS2
    bismark_methylation_extractor --parallel {threads} --no_header --gzip --bedGraph ${{SUBDIR}}/$(basename $SUBDIR)_pe.bam
    """

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class BismarkPEBash(jobs.BashArrayJobMixin, BismarkPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class BismarkMultilanePEBash(jobs.BashArrayJobMixin, BismarkPEBase, jobs.PEFastqBartsMultiLaneMixin):
    pass



class BismarkPEApocrita(BismarkSgeRequirements, BismarkPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class BismarkMultilanePEApocrita(BismarkSgeRequirements, BismarkPEBase, jobs.PEFastqFileIteratorMixin):
    pass


def bismark_run(run_type):
    import argparse
    import sys
    import os

    run_type_dict = {
        'pe_bash': BismarkPEBash,
        'pe_multilane_bash': BismarkMultilanePEBash,
        'pe_apocrita': BismarkPEApocrita,
        'pe_multilane_apocrita': BismarkMultilanePEApocrita,
    }

    if run_type not in run_type_dict:
        raise NotImplementedError("Unrecognised run_type option: %s" % run_type)
    else:
        cls = run_type_dict[run_type]

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-p", "--threads", help="Number of threads", default='1')
    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    required.add_argument("-i", "--index_dir", help="Directory of pre-computed Bismark index", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.read_dir, 'star_alignment')
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
