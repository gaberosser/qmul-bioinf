#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')

from utils import sge


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    optional.add_argument("-o", "--out_dir", help="Output directory", default=os.path.abspath('./'))
    required.add_argument("--project_id", help="SRA project ID (SRX or PRJNA)", required=True)

    # all extra args got to extra
    args, extra = parser.parse_known_args()

    obj = sge.SraDownloadFastqSgeJob(extra_args=extra, **args.__dict__)
    obj.write_submission_script()
    obj.submit()