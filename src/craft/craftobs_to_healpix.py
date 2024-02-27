#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from elasticsearch import Elasticsearch
import healpy as hp
from IPython import embed


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def paginate(e, **kwargs):
    s = e.search(**kwargs)
    while len(s['hits']['hits']) > 0:
        yield s
        # Paginated data from elasticsearch See https://www.elastic.co/guide/en/elasticsearch/reference/8.12/paginate-search-results.html#search-after
        last_sort = s['hits']['hits'][-1]['sort']
        s = e.search(search_after=last_sort, **kwargs)
        
    

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-O','--outfile', help='Output fits file', default='craftscans.fits')
    parser.add_argument('--nside', help='NSIDE for healpix map', type=int, default=200)
    parser.add_argument('--fov-diam', help='FoV diameter to assume. Degrees.', type=float, default=5.5)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    NSIDE = values.nside
    fov_radius = np.radians(values.fov_diam/2) # radians
    resol= hp.nside2resol(NSIDE, arcmin=True)
    npix = hp.nside2npix(NSIDE)
    m = np.zeros(npix, dtype=float)
    print(f'Healpix NSIDE {NSIDE} resolution={resol} arcmin npix={npix}')
    e = Elasticsearch(os.environ['CRAFT_ELASTIC_URL'], ca_certs=os.path.expanduser('~/.config/craco_elastic_ca.crt'))
    index = 'craftscans'
    pit = e.open_point_in_time(index=index, keep_alive='10m')
    pitid = pit['id']
    pitval = {'id':pitid, 'keep_alive':'10m'}

    nscans = 0
    total_duration = 0.0
    try:
        for data in paginate(e, pit=pitval, size=5000, sort='scan_start', source_includes=['ant_ra','ant_dec','duration_hrs']):
            for hit in data['hits']['hits']:
                src = hit['_source']
                theta = 90 - src['ant_dec']
                phi = src['ant_ra']
                theta = src['ant_ra']
                phi = src['ant_dec']
                vec = hp.ang2vec(theta, phi, lonlat=True)
                print(src, theta, phi, vec)
                ipix_disk = hp.query_disc(nside=NSIDE, vec=vec, radius=fov_radius)
                duration = src['duration_hrs']
                m[ipix_disk] += duration
                nscans += 1
                total_duration += duration
    except:
        logging.exception('Exception retrieving data')

    outfile = values.outfile
    title = f'Total duration:{total_duration/24:0.1f} days. {nscans} scans'
    print(title)
    print(f'Writing to {outfile}')
    hp.write_map(outfile, m, overwrite=True, nest=None, coord='C', fits_IDL=False)
    np.save(outfile, m)

    hp.mollview(m, title=title, unit='hrs', coord='C', norm='hist')
    hp.graticule()
    pylab.show()


if __name__ == '__main__':
    _main()
