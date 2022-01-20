#!/usr/bin/env python
"""
Summarises scans in json format to submit to elasticsearch

Copyright (C) CSIRO 2015
"""
import os
import sys
import logging
from .crafthdr import DadaHeader
import glob
import json
from astropy.time import Time
from .sigproc import SigprocFile

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

metafile_cards = ('tsamp','az_start','duration_secs','za_start','nifs','nbits','source','foff','tstart','beamra','beamdec','nchans','data_type','beamno','field_name','fch1','nsamples')

def write_filterbank_meta_file(f, summary, fb, values):
    filname = fb['abspath']
    filshort = os.path.basename(filname)
    metaname = filname+'.meta'
    namespace = 'SB{sbid}/{scanname}/{antname}'.format(**summary)
    asset = {}
    #asset['namespace'] = namespace
    #asset['name'] = filshort
    asset['geoshape/point/latitude'] = fb['beamdec']
    asset['geoshape/point/longitude'] = fb['beamra']
    asset['geoshape/point/elevation'] = 0
    #asset['geoshape/datum'] = 'J2000'

    metafile = open(metaname, 'w+')
    metafile.write('[asset]\n')
    for k, v in asset.items():
        metafile.write('{} = "{}"\n'.format(k, str(v)))


    filterbank = {}
    for c in metafile_cards:
        filterbank[c] = fb[c]

    metafile.write('\n\n[filterbank]\n')
    for k, v in filterbank.items():
        metafile.write('{} = "{}"\n'.format(k, str(v)))
    
    logging.debug('wrote metafile %s', metaname)
    metafile.close()
    

def write_meta_files(f, summary, values):
    for fb in summary['filterbanks']:
        write_filterbank_meta_file(f, summary, fb, values)
    
def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-f','--fname', help='Candidate filename')
    parser.add_argument('-o','--outfile',help='Output json file')
    parser.add_argument('-m','--write-meta-files',help='Write meta files for mediaflux', action='store_true')
    parser.add_argument('-e','--elasticurl', help='Index to this url')

 
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

        if values.elasticurl is not None:
            index_summary(id, summary, values.elasticurl)

        if values.write_meta_files:
            write_meta_files(f, summary, values)
        

    fout.close()

def index_summary(_id, summary, url='akingest01'):
    from elasticsearch import Elasticsearch
    es = Elasticsearch(url)
    es.index(index='craftscans', doc_type='antscan', id=_id, body=summary)
    

def summarise_scan(f):
    hdr = DadaHeader.fromfile(f)
    d = {}
    d['dada_header'] = str(d)
    for key, (v, comment) in hdr.items():
        if 'CHANMAP' in key or key == 'BEAM_ID':
            d[key.lower()] = list(map(int, v.split(',')))
        elif (key.startswith('BEAMFORMER') and 'FREQMAP' in key) or key == "BEAM_RA" or key=='BEAM_DEC':
            d[key.lower()] = list(map(float, v.split(',')))
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
    # convention for elasticsearch is lat, lon, latitude first.
    d['beam_directions'] = list(zip(list(map(fixlon, beam_ra)), beam_dec))
    d['scanname'] = os.path.abspath(f).split('/')[-4]
    assert d['scanname'].startswith('201')

    filfiles = glob.glob(os.path.join(os.path.dirname(f), '*.fil'))
    d['num_filfiles'] = len(filfiles)
    fdetails = []
    for iff, ff in enumerate(filfiles):
        f = {}
        beamno = int(ff.split('.')[-2])
        f['beamno'] = beamno
        f['beamra'] = d['beam_directions'][beamno][0]
        f['beamdec'] = d['beam_directions'][beamno][1]
        f['beam_direction'] = d['beam_directions'][beamno]
        f['abspath'] = os.path.abspath(ff)
        f['atime'] = int(os.path.getatime(ff))
        f['mtime'] = int(os.path.getmtime(ff))
        f['ctime'] = int(os.path.getctime(ff))
        f['filesize'] = int(os.path.getsize(ff))

        fpfile = SigprocFile(ff)
        for hname, value in fpfile.header.items():
            f[hname.lower().replace('.','_')] = value
            
        
        f['data_size_bytes']  = fpfile.data_size_bytes
        f['nsamples'] = fpfile.nsamples
        f['duration_secs'] = fpfile.observation_duration
        f['duration_hrs'] = fpfile.observation_duration/3600.
        f['duration_days'] = fpfile.observation_duration/86400.
        f['source'] = d['source']
        f['field_name'] = d['field_name']
        fdetails.append(f)
            
    d['filterbanks'] = fdetails
    if len(fdetails) > 0:
        tstart = Time(float(fdetails[0]['tstart']), format='mjd', scale='utc')
        d['scan_start_mjd'] = tstart.mjd
        d['scan_start_jd'] = tstart.jd
        d['scan_start'] = tstart.isot
        d['tsamp'] = fdetails[0]['tsamp']
        d['duration_secs'] = fdetails[0]['duration_secs']
        d['duration_hrs'] = fdetails[0]['duration_hrs']
        d['duration_days'] = fdetails[0]['duration_days']
        d['nbits'] =  int(fdetails[0]['nbits'])
        d['nifs'] = int(fdetails[0]['nifs'])
        d['total_filesize'] = sum([f['filesize'] for f in fdetails])

    return d

if __name__ == '__main__':
    _main()
