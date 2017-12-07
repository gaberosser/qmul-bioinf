#!/usr/bin/env python
import os
import argparse
import sys
import math

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')

from pyscript import jobs, sge


class StarSgeRequirements(sge.ApocritaJobMixin):
    @property
    def ram_per_core(self):
        nthread = len(self.params)
        # aim for minimum 32Gb
        gb_per_core = int(math.ceil(32. / float(nthread)))
        return "%dG" % gb_per_core

    @property
    def runtime_mins(self):
        nthread = len(self.params)
        # roughly 60 mins with 12 cores
        return 60 * 12 / float(nthread)


class StarMultilanePEJob(jobs.ArrayJob, jobs.PEFastqIlluminaMultiLaneMixin, StarSgeRequirements):
    title = 'star_multilane_pe'
    required_args = ['read_dir', 'threads', 'genomeDir']
    parameters = [
        # format: (name as it appears in bash script, bash check or None)
        ('$ID', '! -z'),
        ('$READS1', '! -z'),
        ('$READS2', '! -z'),
        ('$SUBDIR', None),  # this will be ignored
        ('$GZ', None)
    ]

    core_cmd = "STAR {extra} --genomeDir {genomeDir} "
    "--readFilesIn $READS1 $READS2 --outFileNamePrefix $SUBDIR --runThreadN {threads}"

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
                ['$(if [ "$GZ" = "gz" ]; then echo "--readFilesCommand zcat" fi)']
            )

    def prepare_submission(self, *args, **kwargs):
        self.setup_params(self.args['read_dir'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-p", "--threads", help="Number of threads", default='1')

    required.add_argument("-i", "--genomeDir", help="Directory of pre-computed STAR index", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.read_dir, 'star_alignment')
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    obj = StarMultilanePEJob(extra_args=extra, **args.__dict__)
    obj.create_script()
    obj.submit()