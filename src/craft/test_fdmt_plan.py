#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from pylab import *
import matplotlib as mpl
from craft import calc11
import numpy as np
from scipy import constants
from craft import uvfits
from craft.craco_plan import PipelinePlan, calc_overlap_channels
from .fdmt_plan import *
from craft.craco import BaselineCell

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def test_fdmt_plan():
    fin = 'b00.uvfits'
    f = uvfits.open(fin)
    plan = PipelinePlan(f, '--ndm 40')
    
    f2 = uvfits.open(fin, skip_blocks=300)
    plan2 = PipelinePlan(f2, '--ndm 40')

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()

    
