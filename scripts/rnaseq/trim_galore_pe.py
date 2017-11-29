#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')

from utils import local_scripting as ls


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("--read_dir", help="Directory containing reads", default='./')
    optional.add_argument("-o", "--out_dir", help="Output directory")

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    if args.out_dir is None:
        # if no output_dir specified, create one in the reads directory
        args.out_dir = os.path.join(args.read_dir, 'trimgalore')
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
        sys.stderr.write("Output directory not specified, using default: %s\n" % args.out_dir)

    obj = ls.TrimgalorePEJob(extra_args=extra, **args.__dict__)
    obj.write_submission_script()
    obj.submit()