import math
from pyscript import jobs, sge
import os


class StarSgeRequirements(sge.ApocritaArrayJobMixin):
    min_ram_gb = 32.
    estimated_runtime_per_core = 720.
    @property
    def ram_per_core(self):
        nthread = int(self.args['threads'])
        # aim for minimum 32Gb
        gb_per_core = int(math.ceil(self.min_ram_gb / float(nthread)))
        print "Num thread: %d" % nthread
        print "Gb per core: %d" % gb_per_core
        return "%dG" % gb_per_core

    @property
    def runtime_mins(self):
        nthread = int(self.args['threads'])
        # roughly 60 mins with 12 cores
        return self.estimated_runtime_per_core / float(nthread)


class StarBaseSE(jobs.ArrayJob):
    title = 'star_se'
    required_args = ['read_dir', 'threads', 'genomeDir']
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$READS1', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
        ('$GZ', None)
    ]
    param_delim = ':'

    core_cmd = "STAR {extra} --genomeDir {genomeDir} " + \
    "--readFilesIn $READS1 --outFileNamePrefix $SUBDIR --runThreadN {threads}"

    def set_default_extras(self):
        if '--outSAMstrandField' not in self.extra_args:
            self.extra_args.extend(
                ['--outSAMstrandField', 'intronMotif']
            )
        if '--quantMode' not in self.extra_args:
            self.extra_args.extend(
                ['--quantMode', 'GeneCounts']
            )
        if '--outSAMtype' not in self.extra_args:
            self.extra_args.extend(
                ['--outSAMtype', 'BAM SortedByCoordinate']
            )
        if '--readFilesCommand' not in self.extra_args:
            self.extra_args.extend(
                ['$([ "$GZ" = "gz" ] && echo "--readFilesCommand zcat" || echo "")']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


class StarBasePE(StarBaseSE):
    title = 'star_pe'
    required_args = ['read_dir', 'threads', 'genomeDir']
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
        ('$GZ', None)
    ]
    param_delim = ':'

    core_cmd = "STAR {extra} --genomeDir {genomeDir} " + \
    "--readFilesIn $READS1 $READS2 --outFileNamePrefix $SUBDIR --runThreadN {threads}"


class StarMultilanePEApocrita(StarSgeRequirements, StarBasePE, jobs.PEFastqIlluminaMultiLaneMixin):
    title = 'star_multilane_pe'
    pass


class StarMultilanePEBash(jobs.BashArrayJobMixin, StarBasePE, jobs.PEFastqIlluminaMultiLaneMixin):
    title = 'star_multilane_pe'
    pass


class StarPEApocrita(StarSgeRequirements, StarBasePE, jobs.PEFastqFileIteratorMixin):
    pass


class StarPEBash(jobs.BashArrayJobMixin, StarBasePE, jobs.PEFastqFileIteratorMixin):
    pass


class StarSEApocrita(StarSgeRequirements, StarBaseSE, jobs.SEFastqFileIteratorMixin):
    pass


class StarSEBash(jobs.BashArrayJobMixin, StarBaseSE, jobs.SEFastqFileIteratorMixin):
    pass


class StarEncodePEBash(jobs.BashArrayJobMixin, StarBasePE, jobs.PEFastqEncodeMultiLaneMixin):
    title = "star_multilane_encode_pe"
    pass


class StarEncodePEApocrita(StarSgeRequirements, StarBasePE, jobs.PEFastqEncodeMultiLaneMixin):
    title = "star_multilane_encode_pe"
    pass


def star_alignment_run(run_type):
    import argparse
    import sys

    run_type_dict = {
        'se_bash': StarSEBash,
        'pe_bash': StarPEBash,
        'pe_multilane_bash': StarMultilanePEBash,
        'pe_encode_bash': StarEncodePEBash,
        'se_apocrita': StarSEApocrita,
        'pe_apocrita': StarPEApocrita,
        'pe_multilane_apocrita': StarMultilanePEApocrita,
        'pe_encode_apocrita': StarEncodePEApocrita,
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

    required.add_argument("-i", "--genomeDir", help="Directory of pre-computed STAR index", required=True)

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
