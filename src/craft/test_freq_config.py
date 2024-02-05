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
from .freq_config import *
from craco.cardcap import CardcapFile
import glob

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def test_from_constructor():
    nc = 6
    mask = np.zeros(nc, dtype=bool)
    mask[0] = True
    f = FrequencyConfig(800.,1.,nc, chan_mask=mask)
    print(str(f))
    assert f.fch1 == 800.
    assert f.foff == 1.
    assert f.nchan == 6
    assert f.fchtop == 805.
    assert f.fcent == 802.5
    assert f.total_bandwidth == 6.0
    assert f.nvalid_channels == 5
    assert f.nmasked_channels == 1
    assert np.all(mask == f.channel_mask)

def test_from_fits_values():
    nc = 6
    mask = np.zeros(nc, dtype=bool)
    mask[2] = True
    f = FrequencyConfig.from_fits_values(fcent=802.5, foff=1., nchan=nc, chan_mask=mask)
    assert f.fch1 == 800.
    assert f.foff == 1.
    assert f.nchan == 6
    assert f.fchtop == 805.
    assert f.fcent == 802.5
    assert f.total_bandwidth == 6.0
    assert f.nvalid_channels == 5
    assert f.nmasked_channels == 1
    assert np.all(mask == f.channel_mask)


def test_fscrunch():
    nc = 6
    mask = np.zeros(nc, dtype=bool)
    mask[0] = True
    forig = FrequencyConfig(800.,1.,nc, chan_mask=mask)
    fscrunch = 3
    f = forig.fscrunch(fscrunch)
    
    assert f.fch1 == forig.freq_of_chan(1)
    assert f.foff == 1.*fscrunch
    assert f.nchan == 6 // fscrunch
    assert f.fchtop == 805. - 1.
    assert f.fcent == 802.5
    assert f.total_bandwidth == 6.0
    assert f.nvalid_channels == 1
    assert f.nmasked_channels == 1
    assert np.all(f.channel_mask == np.array([True, False]))


@pytest.fixture
def cardcap_files():
    pattern = '/data/seren*/big/craco/SB053585/scans/00/2*/ccap*.fits'
    pattern = '/data/seren*/big/craco/SB053553/scans/00/20231006062841/ccap*.fits'

    files = sorted(glob.glob(pattern))
    #ccaps = [open_source(f) for f in files]

    all_files = [CardcapFile(f) for f in files]
    assert len(all_files) > 0
    return all_files

    

def test_from_card_files(cardcap_files):
    f = FrequencyConfig.from_cardcap_files(cardcap_files)
    all_masks = [not f.card_enabled for f in cardcap_files]
    assert f.fch1 == cardcap_files[0].frequencies.min()
    assert abs(f.foff - 0.16666666) < 0.00001
    assert f.nchan == len(cardcap_files)*4
    assert f.nmasked_channels == sum(all_masks)*4


def test_from_card_freqs(cardcap_files):
    all_freqs = [f.frequencies for f in cardcap_files]
    all_masks = [not f.card_enabled for f in cardcap_files]
    all_masks[0] = True

    f = FrequencyConfig.from_freqs_and_masks(all_freqs, all_masks)

    assert f.fch1 == cardcap_files[0].frequencies.min()
    assert abs(f.foff - 0.16666666) < 0.00001
    assert f.nchan == len(cardcap_files)*4
    assert f.nmasked_channels == sum(all_masks)*4


    



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
