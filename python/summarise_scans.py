#!/usr/bin/env python
"""
Summarises scans in json format to submit to elasticsearch

Copyright (C) CSIRO 2015
"""
import os
import sys
import logging
from crafthdr import DadaHeader
import glob
import json
from astropy.time import Time
from sigproc import SigprocFile

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def convert(s):
    if isinstance(s, float) or isinstance(s, int):
        v = s
    elif s.lower() == 'null' or s.lower() == 'none':
        v = None
    elif s is None:
        v = None
    else:
        try:
            v = int(s)
        except:
            try:
                v = float(s)
            except:
                v = s
    return v

def fixlon(l):
    if l > 180:
        l = 180. - l
    
    return l

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-f','--fname', help='Candidate filename')
    parser.add_argument('-o','--outfile',help='Output json file')
 
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

        
    fout = open(values.outfile, 'w')

    for f in values.files:
        summary = summarise_scan(f)
        id = 'SB{sbid}-{scanname}-{antname}'.format(**summary)
        cmd = {'index':{'_index':'craftscans','_type':'antscan', '_id':id}}
        fout.write(json.dumps(cmd))
        fout.write('\n')
        fout.write(json.dumps(summary))
        fout.write('\n')
        logging.info('Summarised scan %s', f)

    fout.close()

def summarise_scan(f):
    hdr = DadaHeader.fromfile(f)
    d = {}
    for key, (v, comment) in hdr.iteritems():
        if 'CHANMAP' in key or key == 'BEAM_ID':
            d[key.lower()] = map(int, v.split(','))
        elif (key.startswith('BEAMFORMER') and 'FREQMAP' in key) or key == "BEAM_RA" or key=='BEAM_DEC':
            d[key.lower()] = map(float, v.split(','))
        elif key == 'FREQ':
            d['freq'] = float(v)
        elif key == 'BEAM_POL':
            d['beam_pol'] = v.split(',')
        else:
            d[key.lower()] = convert(v)

    ra = d['ra']
    dec = d['dec']
    del d['ra']
    del d['dec']
    d['ant_direction'] = [fixlon(ra), dec]
    beam_ra = d['beam_ra']
    beam_dec = d['beam_dec']
    del d['beam_ra']
    del d['beam_dec']
    d['beam_directions'] = zip(map(fixlon, beam_ra), beam_dec)
    d['scanname'] = f.split('/')[-4]
    assert d['scanname'].startswith('201')

    filfiles = glob.glob(os.path.join(os.path.dirname(f), '*.fil'))
    d['num_filfiles'] = len(filfiles)
    fdetails = []
    for ff in filfiles:
        f = {}
        f['abspath'] = os.path.abspath(ff)
        f['atime'] = int(os.path.getatime(ff))
        f['mtime'] = int(os.path.getmtime(ff))
        f['ctime'] = int(os.path.getctime(ff))
        f['size'] = int(os.path.getsize(ff))
        fpfile = SigprocFile(ff)
        for hname, value in fpfile.header.iteritems():
            f[hname.lower().replace('.','_')] = value
            
        f['data_size_bytes']  = fpfile.data_size_bytes
        f['nsamples'] = fpfile.nsamples
        f['duration_secs'] = fpfile.observation_duration
        f['duration_hrs'] = fpfile.observation_duration/3600.
        fdetails.append(f)
            
    d['filterbanks'] = fdetails
    if len(fdetails) > 0:
        tstart = Time(float(fdetails[0]['tstart']), format='mjd', scale='utc')
        d['scan_start_mjd'] = tstart.mjd
        d['scan_start_jd'] = tstart.jd
        d['scan_start'] = tstart.isot
        d['tsamp'] = fdetails[0]['tsamp']
        d['duration_secs'] = fdetails[0]['duration_secs']
        d['nbits'] =  int(fdetails[0]['nbits'])
        d['nifs'] = int(fdetails[0]['nifs'])

    return d

if __name__ == '__main__':
    _main()
