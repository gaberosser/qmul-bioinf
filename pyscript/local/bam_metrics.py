#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import bam_metrics


if __name__ == "__main__":
    bam_metrics.multilane_run("bash")