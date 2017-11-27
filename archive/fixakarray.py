#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import os
import sys
import logging
import json
from fix_headers import TargetParset
from astropy.coordinates import SkyCoord

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def guess_intent(s):
    sbtemplate = s['sb_template']
    field_name = s['field_name']
    scan_intent = 'UNKNOWN'
    if sbtemplate.lower() == 'flyseye':
        if field_name.startswith('PSR'):
            scan_intent = 'PSR_CHECK'
        else:
            scan_intent = 'FRB_FLYSEYE'
    elif sbtemplate.lower() == 'standard':
        if field_name.startswith('PSR') and 'beam' in field_name:
            scan_intent = 'PSR_CAL'
    else:
        logging.error('Unknown sbtemplate %s for field_name %s', sbtemplate, field_name)

    return scan_intent

name_map = {'16:44:49.281,-45:59:09.5':'J1644-4559'}

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--sblist', default='sblist.txt')
    
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    sbdata  = {}
    for line in open(values.sblist, 'rU'):
        if line[0] == '=' or line.startswith('id'):
            continue
        bits = line.split()
        sbdata[int(bits[0])] = bits
    lastkeys = None

    for ifile, fname in enumerate(values.files):
        print 'fname', fname
        foutname = fname+'v2'
        fin = open(fname, 'rU')
        fout = open(foutname, 'w')
        sbid = int(fname.split('.')[0][2:])
        psetname = os.path.join('sbpars_ak','%04d.parset'%sbid)
        pset = TargetParset(psetname)
        for iline, line in enumerate(fin):
            try:
                s = json.loads(line)
            except:
                logging.exception('Could not load line %d of file %s. Skipping', iline, fname)
                continue

            if 'index' in s.keys():
                indexline = line
            else:
                if 'duration_secs' not in s.keys():
                    print 'Ignoring scan', indexline.strip()
                    continue

                s['duration_days'] = s['duration_secs']/86400.
                s['duration_hrs'] = s['duration_secs']/3600.

                ra, dec = s['ant_direction']
                pos = SkyCoord(ra, dec, unit='deg')
                try:
                    src = pset.get_nearest_source(pos)
                    field_name = src['field_name']
                    for oldname, newname in name_map.iteritems():
                        if oldname in field_name:
                            field_name = field_name.replace(oldname, newname)
                            break

                    s['source']  = src['src']
                    s['target'] = field_name
                    s['field_name'] = field_name
                except ValueError, e:
                    logging.exception('Couldnnt find neares source for scan s%s', s)
                    s['source'] = 'UNKNOWN'
                    s['target'] = 'UNKNOWN'
                    s['field_name'] = 'UNKNOWN'

                if 'footprint_pitch' not in s.keys():
                    s['footprint_pitch'] = 0
                    s['footprint_rotation']  = 0
                    s['footprint_shape'] = 'UNKNOWN'


                sbid = s['sbid']
                sdat = sbdata[sbid]
                s['sb_alias'] = sdat[1]
                s['sb_template'] = sdat[3]
                s['sb_template_version'] = sdat[4]
                s['sb_owner'] = sdat[5]
                s['scan_intent'] = guess_intent(s)
                #k = s.keys()
                #if lastkeys is not None:
                    #print 'keydiff', set(lastkeys) - set(k)
                    #print 'k2', set(k) - set(lastkeys)
                    
                #lastkeys = k
                print indexline.strip()
                fout.write(indexline)
                fout.write(json.dumps(s))
                fout.write('\n')
    

if __name__ == '__main__':
    _main()
