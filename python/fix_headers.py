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
import datetime
import glob
import re
import shutil
from astropy.coordinates import SkyCoord
from aces.footprint_class import Footprint
from aces.skypos import skypos

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def mysep(pos1, pos2):
    ''' Returns separation in radians'''
    sep =np.arccos(np.sin(pos1.dec.rad)*np.sin(pos2.dec.rad) + np.cos(pos1.dec.rad)*np.cos(pos2.dec.rad)*np.cos(pos1.ra.rad - pos2.ra.rad))
    return sep

def mysep_skypos(pos1, pos2):
    sep =np.arccos(np.sin(pos1.dec)*np.sin(pos2.dec) + np.cos(pos1.dec)*np.cos(pos2.dec)*np.cos(pos1.ra - pos2.ra))
    return sep
        
class TargetParset(dict):
    def __init__(self, f):

        self.sources = {}
        self.generic_target_parms = {}
        self._parse_sources(f)
        self._add_field_directions()

    def _parse_sources(self, pset):
        sources = self.sources
        with open(pset, 'rU') as f:
            for iline, line in enumerate(f):
                if line.startswith('#'):
                    continue
                if line.startswith('='):
                    continue
                if not '=' in line:
                    continue

                key, value = line.split('=')
                key = key.strip()
                self[key] = value.strip()
                keyparts = key.strip().split('.')
                if key.startswith('common.target.src') and '%d' not in key:
                    src = keyparts[2]

                    if src not in sources:
                        sources[src] = {'src':src}

                    srckey = keyparts[3]
                    sources[src][srckey] = value.strip()

        return sources
        
    def get_target_param(self, src, param):
        assert '%' not in src
        assert src.startswith('src')
        p = None

        targ_str = 'common.target.{}.{}'.format(src, param)
        if targ_str in self.keys():
            p = self[targ_str]
        else:
            srcid = int(src[3:])
            for k, v in self.iteritems():
                if '%d' not in k:
                    continue

                kreplaced = k % srcid
                if kreplaced == targ_str:
                    p = v
                    break
                
        if p is None:
            raise ValueError('Couldnt find target src {} param {}'.format(src, param))
        
        return p

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
            pastr = self.get_target_param(src, 'pol_axis')
            pol_type, pol_angle = pastr.replace('[','').replace(']','').split(',')
            data['pol_type'] = pol_type
            data['pol_angle'] = float(pol_angle)

            fp = Footprint.named('square_6x6', np.radians(0.9), np.radians(float(pol_angle)+45))
            fp.setRefpos(sp)
            data['footprint'] = fp

    def get_nearest_source(self, pos, tol_deg=0.1):
        sources = [s for s in self.sources.values() if 'skycoord' in s]
        nearest = min(sources, key=lambda p: mysep(p['skycoord'], pos))
        assert mysep(nearest['skycoord'], pos) <= np.radians(tol_deg)

        return nearest

def guess_intent(hdr):
    sbtemplate = hdr.get_value('SB_TEMPLATE')
    field_name = hdr.get_value('FIELD_NAME')
    scan_intent = 'UNKNOWN'
    if sbtemplate.lower() == 'flyseye':
        if field_name.startswith('G'):
            scan_intent = 'FRB_FLYSEYE'
        elif field_name.startswith('PSR'):
            scan_intent = 'PSR_CHECK'
        else:
            logging.error('Could not work out intent for sbtempalte %s field_name %s', sbtemplate, field_name)
    elif sbtemplate.lower() == 'standard':
        if field_name.startswith('PSR') and 'beam' in field_name:
            scan_intent = 'PSR_CAL'
    else:
        logging.error('Unknown sbtemplate %s for file %s', sbtemplate, hdr_file)

    hdr += ('SCAN_INTENT', scan_intent, 'Scan intent')

name_map = {'16:44:49.281,-45:59:09.5':'J1644-4559'}


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Fix headers', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-p', '--parset', help='Parset to get field names from')
    parser.add_argument('-l','--sblist', help='Sschedblock list to get data from', default='sbpars/sblist.txt')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    pset = None
    if values.parset is not None:
        pset = TargetParset(values.parset)
    
    sbdata = {}
        
    for line in open(values.sblist, 'rU'):
        if line[0] == '=' or line.startswith('id'):
            continue
        bits = line.split()
        sbdata[int(bits[0])] = bits


    for hdr_file in values.files:
        try:
            fix_header(hdr_file, sbdata, pset, values)
        except:
            logging.exception('Could not fix header file %s', hdr_file)
            raise

def fix_header(hdr_file, sbdata, pset, values):
    hdr = DadaHeader.fromfile(hdr_file)

    sbid = int(hdr['SBID'][0])
    
    if pset is None:
        pset = TargetParset(os.path.join('sbpars', 'SB%04d.parset'%sbid))

    sbinfo = sbdata[sbid]
    
    sbtemplate = sbinfo[3]
    if 'FIXED_VERSION' in hdr.keys():
        last_version = int(hdr.get_value('FIXED_VERSION'))
        hdr.set_value('FIXED_VERSION', last_version+1)
    else:
        hdr += 'FIXED_VERSION', 1, 'Fix version - increments by 1 on each change'

    hdr += ('FIXED_UTC', datetime.datetime.utcnow().isoformat(), 'Dated fixed')
    hdr += ('SB_ALIAS', sbinfo[1],'SB alias')
    hdr += ('SB_TEMPLATE', sbtemplate, 'SB template')
    hdr += ('SB_TEMPLATE_VERSION', sbinfo[4], 'SB template version')
    hdr += ('SB_OWNER', sbinfo[5], 'SB owner')

    if sbtemplate.lower() != 'beamform':
        ra, dec = map(float, (hdr['RA'][0], hdr['DEC'][0]))
        pos = SkyCoord(ra, dec, unit='deg')
        src = pset.get_nearest_source(pos)
        fp = src['footprint']
        # some headers had incorrect positions - not sure which now. Some headers don't have footprints either
        #hdr.set_value('BEAM_RA', ','.join(map(str, np.degrees([p.ra for p in fp.getPositions()]))))
        #hdr.set_value('BEAM_DEC', ','.join(map(str, np.degrees([p.dec for p in fp.getPositions()]))))

        #hdr.set_value('RA', src['skycoord'].ra.deg)
        #hdr.set_value('DEC',  src['skycoord'].dec.deg)

        for f in ('field_direction',):
            hdr.set_value(f.upper(), src[f])

        field_name = src['field_name']

        # sometimes names were incorrect in headers. Fix using the name map
        for oldname, newname in name_map.iteritems():
            if oldname in field_name:
                field_name = field_name.replace(oldname, newname)
                break
            
        hdr.set_value('SOURCE', src['src'])
        hdr.set_value('FIELD_NAME', field_name)
        hdr.set_value('TARGET', field_name)

        guess_intent(hdr)

    shutil.move(hdr_file, hdr_file+'.beforefix')
    
    with open(hdr_file, 'w') as fout:
        hdr.tofile(fout, add_zeros=False)
        

if __name__ == '__main__':
    _main()
