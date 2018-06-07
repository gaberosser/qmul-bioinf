import math
import os
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
        ('$NAME', None),
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
    ]
    param_delim = ':'
    core_cmd = "salmon quant -i {index_dir} -l {library_type} -p {threads} -1 $READS1 -2 $READS2 -o $SUBDIR {extra}"

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
        ('$NAME', None),
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


class SalmonMultilaneEncodePEApocrita(SalmonSgeRequirements, SalmonPEBase, jobs.PEFastqEncodeMultiLaneMixin):
    file_sep = ' '  # the character used to separate files of the same read number in different lanes
    pass


class SalmonMultilaneEncodePEBash(jobs.BashArrayJobMixin, SalmonPEBase, jobs.PEFastqEncodeMultiLaneMixin):
    file_sep = ' '  # the character used to separate files of the same read number in different lanes
    pass


class SalmonMultilaneIlluminaPEApocrita(SalmonSgeRequirements, SalmonPEBase, jobs.PEFastqIlluminaMultiLaneMixin):
    file_sep = ' '  # the character used to separate files of the same read number in different lanes
    pass


class SalmonMultilaneIlluminaPEBash(jobs.BashArrayJobMixin, SalmonPEBase, jobs.PEFastqIlluminaMultiLaneMixin):
    file_sep = ' '  # the character used to separate files of the same read number in different lanes
    pass



class SalmonPEApocrita(SalmonSgeRequirements, SalmonPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class SalmonPEBash(jobs.BashArrayJobMixin, SalmonPEBase, jobs.PEFastqFileIteratorMixin):
    pass


class SalmonSEApocrita(SalmonSgeRequirements, SalmonSEBase, jobs.SEFastqFileIteratorMixin):
    pass


class SalmonSEBash(jobs.BashArrayJobMixin, SalmonSEBase, jobs.SEFastqFileIteratorMixin):
    pass


def run_salmon(run_type):
    import argparse
    import sys

    run_type_dict = {
        'se_bash': SalmonSEBash,
        'pe_bash': SalmonPEBash,
        'pe_multilane_bash': SalmonMultilaneIlluminaPEBash,
        'pe_multilane_encode_bash': SalmonMultilaneEncodePEBash,
        'se_apocrita': SalmonSEApocrita,
        'pe_apocrita': SalmonPEApocrita,
        'pe_multilane_apocrita': SalmonMultilaneIlluminaPEApocrita,
        'pe_multilane_encode_apocrita': SalmonMultilaneEncodePEApocrita,
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
    optional.add_argument("-l", "--library_type", help="Library type (default: auto)", default='A')
    optional.add_argument("--include", help="List of filestems to include (comma separated)")
    optional.add_argument("--exclude", help="List of filestems to exclude (comma separated)")

    required.add_argument("-i", "--index_dir", help="Directory of pre-computed Salmon index", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.read_dir, 'salmon')
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
