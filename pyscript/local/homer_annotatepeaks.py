#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import homer


if __name__ == "__main__":
    homer.homer_annotate_peaks_run("bash")