#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2022
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import pytest
from craft import craco_plan
from astropy.coordinates import SkyCoord, Angle
from craft.craco_plan import ImageParams1d, ImageParams2d, splitn

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def test_image_params_1d_os():


    # Test sample where you hit the limit in oversampling
    lam = 0.2 #m
    bmax = 2e3/lam # wavelengths
    npix = 256
    os = 2.1
    fov = Angle('1.05deg')
    p1 = ImageParams1d(bmax, npix, os, fov)
    print(p1)
    assert p1.os == os
    assert p1.fov < fov

def test_image_params_1d_fov():
    # Test sample where you achieve the fov
    lam = 0.2 #m
    bmax = 1.3e3/lam # wavelengths
    npix = 256
    os = 2.1
    fov = Angle('1.05deg')
    p1 = ImageParams1d(bmax, npix, os, fov)
    print(p1)
    assert p1.os >= os
    assert p1.fov == fov

def test_image_params_1d_1ghz():
    # Test sample where you achieve the fov
    fmax = 1e9
    c = 3e8
    lam = c/fmax
    bmax = 1.3e3/lam # wavelengths
    npix = 256
    os = 2.1
    fov = Angle('1.05deg')
    p1 = ImageParams1d(bmax, npix, os, fov)
    print(p1)
    assert p1.os >= os
    assert p1.fov == fov

def test_image_params_2d():
    fmax = 1e9
    c = 3e8
    lam = c/fmax
    bmax = 1.3e3/lam # wavelengths
    npix = 256
    os_str = '2.1'
    fov_str = '1.05deg'
    uvmax = (bmax,bmax*0.9)
    p2 = ImageParams2d(uvmax, npix, os_str, fov_str)
    print(p2)
    print('LMCELL', p2.lmcell, type(p2.lmcell))
    print('UVCELL', p2.uvcell, type(p2.uvcell))


def test_splitn():
    assert splitn('a',',',2) == ['a','a']
    assert splitn('a,b',',',2) == ['a','b']

    

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
