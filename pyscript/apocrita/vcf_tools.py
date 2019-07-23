#!/usr/bin/env python
import os
import argparse
import sys
import math

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import vcf_tools


if __name__ == "__main__":
    vcf_tools.run("apocrita")