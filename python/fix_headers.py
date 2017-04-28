#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
from crafthdr import DadaHeader
import glob
import re
from astropy.coordinates import SkyCoord
from aces.footprint_class import Footprint
from aces.skypos import skypos

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def mysep(pos1, pos2):
    sep =np.arccos(np.sin(pos1.dec.rad)*np.sin(pos2.dec.rad) + np.cos(pos1.dec.rad)*np.cos(pos2.dec.rad)*np.cos(pos1.ra.rad - pos2.ra.rad))
    return sep

        
class TargetParset(object):
    def __init__(self, f):

        self.sources = {}
        self._parse_sources(f)
        self._add_field_directions()

    def _parse_sources(self, pset):
        sources = self.sources
        with open(pset, 'rU') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if not '=' in line:
                    continue

                key, value = line.split('=')
                key = key.strip()
                keyparts = key.strip().split('.')
                if key.startswith('common.target.src'):
                    src = keyparts[2]

                    if src not in sources:
                        sources[src] = {'src':src}

                    srckey = keyparts[3]
                    sources[src][srckey] = value.strip()

        return sources

    def _add_field_directions(self):
        for src, data in self.sources.iteritems():
            if 'field_direction' not in data:
                continue

            fdstr = data['field_direction']
            fdbits = fdstr.replace('[','').replace(']','').split(',')
            rastr, decstr, epochstr = fdbits
            skycoord = SkyCoord(rastr, decstr, unit='hour, deg')
            sp = skypos(rastr, decstr)
            data['skycoord'] = skycoord
            data['skypos'] = sp

            print fdstr, rastr, decstr, skycoord, sp

            pastr = data['pol_axis']
            pol_type, pol_angle = pastr.replace('[','').replace(']','').split(',')
            data['pol_type'] = pol_type
            data['pol_angle'] = float(pol_angle)

            fp = Footprint.named('square_6x6', np.radians(0.9), np.radians(float(pol_angle)+45))
            fp.setRefpos(sp)
            data['footprint'] = fp

    def get_nearest_source(self, pos):
        sources = [s for s in self.sources.values() if 'skycoord' in s]
        nearest = min(sources, key=lambda p: mysep(p['skycoord'], pos))
        #for p in sources:

        # GOLLY, POL.SEPARATION IS *SLOW*
        #print p['field_name'], p['skycoord'], pos, pos.separation(p['skycoord']), np.degrees(mysep(pos, p['skycoord']))
        #sources = [s for s in self.sources.values() if 'skypos' in s]

        #str(pos) # need this to populate pos.ras and pos.decs

        #for p in sources:
        #    print pos, p['field_name'], p['skypos'], p['skypos'].distDeg(pos.ra, pos.dec)

        #nearest = min(sources, key=lambda p: p['skypos'].dist(pos.ra, pos.dec))

        return nearest

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-p', '--parset', help='Parset to get field names from')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    pset = TargetParset(values.parset)

    for hdr_file in values.files:
        hdr = DadaHeader.fromfile(hdr_file)
        ra, dec = map(float, (hdr['RA'][0], hdr['DEC'][0]))
        pos = SkyCoord(ra, dec, unit='deg')
        src = pset.get_nearest_source(pos)
        hdr.set_value('RA', src['skycoord'].ra.deg)
        hdr.set_value('DEC',  src['skycoord'].dec.deg)
        fp = src['footprint']
        hdr.set_value('BEAM_RA', ','.join(map(str, np.degrees([p.ra for p in fp.getPositions()]))))
        hdr.set_value('BEAM_DEC', ','.join(map(str, np.degrees([p.dec for p in fp.getPositions()]))))
        
        for f in ('field_direction',):
            hdr.set_value(f.upper(), src[f])

        hdr.set_value('SOURCE', src['field_name'])

        with open(hdr_file+'.fixed', 'w') as fout:
            hdr.tofile(fout, add_zeros=False)
        

if __name__ == '__main__':
    _main()
