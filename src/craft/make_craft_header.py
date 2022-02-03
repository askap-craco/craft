#!/usr/bin/env python
"""
Sets up craft and downlaods channel mapping

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys
import askap.craft.beamformer 
from askap.craft.beamformer import CraftBeamformer, FilterTypes, ZoomModes
from askap.craft.freqconfig import FreqConfig
import numpy as np
from multiprocessing import Pool

def _main():
    from argparse import ArgumentParser
    from askap.craft.cmdline import strrange
    parser = ArgumentParser(description='Sets up craft and downloads channel mapping')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-a', '--antennas', help='Antenna number, comma sparated or dash. e.g. "3", "13-15", "1,2,14"', type=strrange)
    parser.add_argument('-c','--cards', help='Card number. comma separated, or dash. e.g. "3", "1-7" or "1,2,4-7"', default="1-7", type=strrange)
    parser.add_argument('-i','--int-time', help='Integration time (1-65535)', default=1000, type=int)
    parser.add_argument('-s','--int-cycles', help='Number of cycles to combine (1-7)', default=7, type=int)
    parser.add_argument('-n','--beam-number', help='Beam number to save voltages for (0-35). Unspecified for all beams', default=None)
    parser.add_argument('-m','--buffer-mode', help='Buffer number of bits', choices=[16, 8, 4, 1], type=int, default=16)
    parser.add_argument('-b','--band',help='Band number 0 = 1200 MHz BPF, 1 = 1450 MHz BPF, 2 - 1800 MHz BPF, 3 - LPF', choices=[0,1,2,3], type=int)
    parser.add_argument('-z','--zoom',help='Zoom mode', choices=[0,1,2,3,4,5], default=0, type=int)
    parser.add_argument('-f','--center-freq', help='Center frequency (MHz)', type=int)
    parser.add_argument('-d','--pushd-delay', help='PUsh download delay -1=auto, 0=none, x=value (seconds)', type=float, default=0)
    parser.add_argument('--program-fpga', help="Force programming of FPGA. 0=Don't. 1=Do", type=int, choices=[0,1], default=0)
    parser.add_argument('-o','--freqfile', help='File to write frequency maping to', default='.freqs')
    parser.add_argument('--corrblock', help='Correlator block to which to apply center freq, unspecified or 0 for default', type=int, choices=[0,1,2,3,4,5,6,7,8], default=0)
    parser.add_argument('--num-threads', help='Number of threads to spawn', type=int, default=1)
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    cards = values.cards
    ants = values.antennas
    num_beamformers = len(cards)*len(ants)

    print('NUM_BEAMFORMERS {}'.format(num_beamformers))
    for iant, ant in enumerate(ants):
        for icard, c in enumerate(cards):
            bfno = icard + iant*len(cards)
            print('BEAMFORMER{}_ADDR 10.2.{}.{}'.format(bfno, ant, c))
            chanmap = ','.join(map(str, np.arange(48)))
            print('BEAMFORMER{}_CHANMAP {}'.format(bfno, chanmap))
            print('BEAMFORMER{}_FREQMAP {}'.format(bfno, chanmap))
            

            

if __name__ == '__main__':
    _main()
