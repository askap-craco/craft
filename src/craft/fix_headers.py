#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
from .crafthdr import DadaHeader
import datetime
import glob
import re
import shutil
from astropy.coordinates import SkyCoord
from aces.footprint_class import Footprint
from aces.skypos import skypos

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def mysep2(pos1, pos2):
    return pos1.separation(pos2).rad
        
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
        if targ_str in list(self.keys()):
            p = self[targ_str]
        else:
            srcid = int(src[3:])
            for k, v in self.items():
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
        for src, data in self.sources.items():
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

            #fp = Footprint.named('square_6x6', np.radians(0.9), np.radians(float(pol_angle)+45))
            #fp.setRefpos(sp)
            #data['footprint'] = fp

    def get_nearest_source(self, pos, tol_deg=0.1):
        sources = [s for s in list(self.sources.values()) if 'skycoord' in s]
        nearest = min(sources, key=lambda p: mysep2(p['skycoord'], pos))
        sepdeg  = np.degrees(mysep2(nearest['skycoord'], pos) )
        if sepdeg > tol_deg:
            raise ValueError('Nearest source to pos {} was {} but separation was {}'.format(pos, nearest, sepdeg))

        return nearest

def guess_intent(hdr):
    sbtemplate = hdr.get_value('SB_TEMPLATE')
    field_name = hdr.get_value('FIELD_NAME')
    scan_intent = 'UNKNOWN'
    if 'flyseye' in sbtemplate.lower():
        if field_name.startswith('PSR'):
            scan_intent = 'PSR_CHECK'
        else:
            scan_intent = 'FRB_FLYSEYE'
    elif sbtemplate.lower() == 'standard':
        if field_name.startswith('PSR') and 'beam' in field_name:
            scan_intent = 'PSR_CAL'
    else:
        logging.error('Unknown sbtemplate %s for ', sbtemplate)
                      
    hdr += ('SCAN_INTENT', scan_intent, 'Scan intent')

name_map = {'16:44:49.281,-45:59:09.5':'J1644-4559'}


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Fix headers', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-p', '--parset', help='Parset to get field names from')
    parser.add_argument('-l','--sblist', help='Sschedblock list to get data from')
    parser.add_argument('-d','--parset-dir', help='Directory containing parset', default='sbpars')
    parser.add_argument('-t','--tolerance', help='search tolerance when associating a header RA/DEC with a parset RA/Dec (degrees)', type=float, default=0.1)
    parser.add_argument('--fix', action='store_true', help='Actually fix the header - otherwise dont change anything')
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

    if values.sblist is not None:
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
    logging.info('Fixing header %s', hdr_file)

    sbid = int(hdr['SBID'][0])
    
    if pset is None:
        pset = TargetParset(os.path.join(values.parset_dir, 'SB%05d.parset'%sbid))

    if len(sbdata) > 0: # if sbinfo is there, we should use it and fail if the SBID isn't in there.
        sbinfo = sbdata[sbid]
    else:
        sbinfo = None # if there's no sbdata at all, then don't complain
    
    if 'FIXED_VERSION' in list(hdr.keys()):
        last_version = int(hdr.get_value('FIXED_VERSION'))
        new_version = last_version+1
        hdr.set_value('FIXED_VERSION', new_version)
    else:
        new_version = 2
        hdr += 'FIXED_VERSION', new_version, 'Fix version - increments by 1 on each change'


    hdr += ('FIXED_UTC', datetime.datetime.utcnow().isoformat(), 'Dated fixed')
    if sbinfo is None:
        sbtemplate = hdr['SB_TEMPLATE'][0]
    else:
        sbtemplate = sbinfo[3]
        hdr += ('SB_ALIAS', sbinfo[1],'SB alias')
        hdr += ('SB_TEMPLATE', sbtemplate, 'SB template')
        hdr += ('SB_TEMPLATE_VERSION', sbinfo[4], 'SB template version')
        hdr += ('SB_OWNER', sbinfo[5], 'SB owner')

    if sbtemplate.lower() != 'beamform':
        ra, dec = list(map(float, (hdr['RA'][0], hdr['DEC'][0])))
        ant_pos = SkyCoord(ra, dec, unit='deg')
        src = pset.get_nearest_source(ant_pos, values.tolerance)
        field_direction = src['field_direction']

        field_ra, field_dec, field_epoch = field_direction.replace("'",'').replace('[','').replace(']','').replace(' ','').split(',')
        
        field_pos = SkyCoord(field_ra + ' ' + field_dec, unit=('hourangle','deg'))
        logging.debug('Nearest src %s', src)
        logging.debug('Field pos %s = ra %s dec %s, decimal=%s', field_pos.to_string('hmsdms'), field_ra, field_dec, field_pos.to_string('decimal'))
        logging.debug('Ant pos %s %s Separtion to requrested pos = %sarcsec', ant_pos.to_string('hmsdms'), ant_pos.to_string('decimal'), field_pos.separation(ant_pos).arcsec)
        assert field_pos.separation(ant_pos).arcsec < 100

        for f in ('field_direction',):
            hdr.set_value(f.upper(), src[f])

        field_name = src['field_name']

        # sometimes names were incorrect in headers. Fix using the name map
        for oldname, newname in name_map.items():
            if oldname in field_name:
                field_name = field_name.replace(oldname, newname)
                break
            
        hdr.set_value('SOURCE', src['src'])
        hdr.set_value('FIELD_NAME', field_name)
        hdr.set_value('TARGET', field_name)
        hdr += 'TARGET_POL_ANGLE', src['pol_angle'], 'Target polarisation angle'
        hdr += 'TARGET_POL_TYPE', src['pol_type'], 'Target polarisation tracking'

        ant_pa = float(hdr['TARGET_POL_ANGLE'][0])
        antpos = SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg')
        fp_shape = hdr['FOOTPRINT_NAME'][0]
        fp_pitch = float(hdr['FOOTPRINT_PITCH'][0])
        fp_pa = float(hdr['FOOTPRINT_ROTATION'][0])
        fp = Footprint.named(fp_shape, np.radians(fp_pitch), np.radians(fp_pa + ant_pa))
        fp.setRefpos(np.radians([field_pos.ra.deg, field_pos.dec.deg]))
        logging.debug('FOOTPRINT %s %s %s %s %s %s %s %s', fp_shape, fp_pitch, fp_pa, 'ant pa', ant_pa, 'refpos', ra, dec)
        fp_ras = np.degrees([p.ra for p in fp.positions])
        fp_decs = np.degrees([p.dec for p in fp.positions])
        ra_strings = list(map(str, fp_ras))
        ra_strings += ra_strings
        dec_strings = list(map(str, fp_decs))
        dec_strings += dec_strings
        hdr.set_value('BEAM_RA',','.join(ra_strings))
        hdr.set_value('BEAM_DEC',','.join(dec_strings))
        hdr.set_value('RA', field_pos.ra.deg)
        hdr.set_value('DEC', field_pos.dec.deg)
        hdr += 'ANT_RA', ant_pos.ra.deg, 'Antenna actual RA at the beginning of the scan (unreliable by a few arcmin)'
        hdr += 'ANT_DEC', ant_pos.dec.deg, 'Antenna actual dec at the beginning of the scan (unreliable by a few arcmin)'
        

        
        
        

        guess_intent(hdr)

    #shutil.move(hdr_file, hdr_file+'.beforefix')
    outfile = '{}.v{}'.format(hdr_file, new_version)
    
    if values.fix:
        with open(outfile, 'w') as fout:
            hdr.add_comment = False
            logging.info('Writing new header to %s', outfile)
            hdr.tofile(fout, add_zeros=False)
        

if __name__ == '__main__':
    _main()
