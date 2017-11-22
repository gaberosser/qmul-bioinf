#!/usr/bin/env python
import os
import re
import subprocess
import argparse
import sys
import datetime
import csv

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')

from utils.log import get_file_logger
from utils.output import unique_output_dir

TITLE = 'cufflinks'
CORE_CMD = 'cufflinks'
LOG_DIR = os.path.join(os.environ['HOME'], 'log')

# list of regexes to clean up output dir. Each entry gives (pattern, replace) inputs to re.sub
# NB we will already have stripped the .bam file extension
outdir_cleanup_regexes = [
    (r'Aligned\.sortedByCoord\.out', ''),
    (r'Aligned\.out', ''),
]

if __name__ == "__main__":
    """
    Usage: cufflinks.py path_to_bams path_to_gtf path_to_output_dir passed_to_function
    """
    now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    logger = get_file_logger(
        name=__name__,
        filestem=os.path.join(
            os.environ['HOME'],
            'log',
            'cufflinks.%s.log' % now_str
        )
    )

    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")
    optional.add_argument("-p", "--threads", help="Number of threads", default='1')
    optional.add_argument("--library_type", help="Library type", default='fr-unstranded')

    required.add_argument("-G", "--GTF", help="GTF annotations for quantification", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        test_out = os.path.join(args.read_dir, 'cufflinks')
        if not os.path.isdir(test_out):
            os.makedirs(test_out)
            args.out_dir = test_out
        else:
            args.out_dir = unique_output_dir('cufflinks', root_output_dir=args.read_dir)
        sys.stderr.write("Created output directory %s\n" % args.out_dir)

    read_dir = args.read_dir
    ref_fn = args.GTF
    out_dir = args.out_dir
    threads = args.threads
    library_type = args.library_type

    # write arguments to a file for reference
    conf_fn = os.path.join(out_dir, "%s.conf" % now_str)
    with open(conf_fn, 'w') as f:
        c = csv.writer(f, delimiter='\t')
        c.writerow(['Field', 'Value'])
        c.writerows([(k, str(v)) for k, v in args.__dict__.items()])

    logger.info(
        "Initialised logger from directory %s. However, the working directory for cluster work is %s.",
        os.path.abspath(os.path.curdir),
        out_dir,
    )

    logger.info("Arguments: read_dir=%s, ref_fn=%s, out_dir=%s, threads=%s, library_type=%s",
                read_dir, ref_fn, out_dir, threads, library_type)

    if not os.path.isfile(ref_fn):
        logger.error("Reference file %s does not exist", ref_fn)
        raise AttributeError("Unable to find reference file.")

    rr = re.compile(r'\.bam$', flags=re.IGNORECASE)
    flist = [t for t in os.listdir(read_dir) if re.search(rr, t)]
    # check for existing output and identify pairs of files
    fl = {}
    for t in flist:
        base = re.sub(r'\.bam', '', t)
        # apply cleanup
        for patt, repl in outdir_cleanup_regexes:
            base = re.sub(patt, repl, base)
        out_subdir = os.path.abspath(os.path.join(out_dir, base))
        logger.info("Input file %s. Cleaned filestem %s. Output subdir %s.", t, base, out_subdir)
        # if SAM output file exists, log warning and skip
        if os.path.isdir(out_subdir):
            if len(os.listdir(out_subdir)) > 0:
                logger.warn("Dir already exists: %s. Skipping.", out_subdir)
                continue
            else:
                logger.info("Using existing empty output subdir %s", out_subdir)
        else:
            logger.info("Created output subdir %s", out_subdir)
            os.makedirs(out_subdir)

        fl[base] = {
            'bam': os.path.abspath(os.path.join(read_dir, t)),
            'outdir': out_subdir,
        }

    # log the filelist
    logger.info("Found %d BAM files: %s.", len(fl), ', '.join(fl.keys()))

    # create a text file containing the arguments that will allow us to request an array
    param_fn = os.path.join(out_dir, "params.%s.txt" % now_str)
    with open(param_fn, 'wb') as f:
        f.writelines(["%s,%s\n" % (t['bam'], t['outdir']) for t in fl.values()])

    # generate the command to run
    cmd = "cufflinks --no-update-check -G {ref_fn} -p {threads} -o $OUTDIR --library-type {libtype} $BAM".format(
        ref_fn=ref_fn,
        threads=threads,
        libtype=library_type
    )

    # create the actual submission script
    script_fn = os.path.join(out_dir, "cufflinks.%s.sh" % now_str)
    with open(script_fn, 'wb') as f:
        sh = \
"""
#!/bin/sh
#$ -wd {work_dir} # Set the working directory for the job to the current directory
#$ -j y           # Join stdout and stderr
#$ -pe smp {threads}      # Request {threads} CPU cores
#$ -l h_rt=1:0:0  # Request 1 hour runtime
#$ -l h_vmem=512M   # Request .5GB RAM / core
#$ -t 1-{nfile}
# source bashrc to get correct path
. $HOME/.bashrc
# define params
INPUTS=$(sed -n "${{SGE_TASK_ID}}p" {params_fn})
BAM=$(echo $INPUTS | cut -d , -f 1)
OUTDIR=$(echo $INPUTS | cut -d , -f 2)

# actual code
if [[ -f $BAM && ! -z $OUTDIR ]]; then
{cmd}
else
echo "Unable to execute run ${{SGE_TASK_ID}} as the read file did not exist or the output dir variable is empty."
echo "Read file: $BAM"
echo "Output dir: $OUTDIR"
fi
""".format(threads=threads, nfile=len(fl), params_fn=param_fn, cmd=cmd, work_dir=out_dir)
        f.write(sh)

    logger.info("Cluster submission script written: %s", sh)
    subprocess.call(['qsub', script_fn])
