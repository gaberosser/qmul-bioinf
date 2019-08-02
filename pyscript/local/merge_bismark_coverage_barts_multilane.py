#!/usr/bin/env python
import os
import argparse
import sys
import math

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import merge_bismark_coverage


if __name__ == "__main__":
    merge_bismark_coverage.run("bash_multilane_barts")