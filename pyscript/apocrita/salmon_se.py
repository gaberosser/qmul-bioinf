#!/usr/bin/env python
import os
import argparse
import sys
import math

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import salmon


if __name__ == "__main__":
    salmon.run_salmon("se_apocrita")