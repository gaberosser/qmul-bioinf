#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import bowtie2


if __name__ == "__main__":
    bowtie2.multilane_run("bash")