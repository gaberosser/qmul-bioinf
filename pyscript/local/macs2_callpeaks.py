#!/usr/bin/env python
import os
import argparse
import sys

# add root of project dir to the path
sys.path.append(os.path.dirname(__file__) + '/../../')
from pyscript import macs2


if __name__ == "__main__":
    macs2.callpeaks_multilane_run("bash")