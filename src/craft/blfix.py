#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
from .cmdline import strrange
from .craftcor import CLIGHT

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-a','--antennas', type=strrange, help='ASKAP antenna numbers')
    parser.add_argument('--fcm', help='Original FCM')
    parser.add_argument('--blfile', help='blfit output')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    pset = {}
    for line in open(values.fcm,'rU'):
        line = line.strip()
        if line.startswith('#'):
            continue
        bits = line.split('=')
        if len(bits) == 2:
            pset[bits[0].strip()] = bits[1].strip()

    blin = open(values.blfile, 'rU').read()
    blidx = blin.find('Miriad baseline convention')
    print('# blfix.py args', ' '.join(sys.argv))
    
    for line in blin[blidx:].split('\n'):
        if 'Miriad' in line or 'Ant' in line or '----' in line:
            continue

        line = line.strip()
        if line == '':
            break

        iant = int(line[0]) # miriad 1-index antenna number
        restline = line[1:]        # Fortran fixed widths aren't big enough
        dxyz = np.array(list(map(float, restline.split())))
        dxyz_meters = dxyz*CLIGHT/1e9
        antno = values.antennas[iant-1]
        s = 'common.antenna.ant{:d}.location.itrf'.format(antno)
        oldxyz = pset[s]
        oldxyz = np.array(list(map(float, oldxyz[1:-1].split(','))))
        newxyz = oldxyz + dxyz_meters
        sout = '{} = [{:0.4f}, {:0.4f}, {:0.4f}]'.format(s, newxyz[0],newxyz[1],newxyz[2])
        print('# iant={} antno={} dXYZ={} m oldxyz={} m'.format(iant, antno, dxyz_meters, oldxyz))
        print(sout)
        
        
        
        

        

    

    
    

if __name__ == '__main__':
    _main()
