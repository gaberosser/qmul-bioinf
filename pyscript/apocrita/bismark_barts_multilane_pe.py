#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import bismark


if __name__ == "__main__":
    bismark.bismark_run('pe_multilane_apocrita')