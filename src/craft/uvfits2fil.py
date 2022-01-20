#!/usr/bin/env python
"""
conver uvfits file to filterbank
Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from . import sigproc
from astropy.io import fits

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def get_antname( hdu, iant):
    an = hdu['AIPS AN']
    return an.data['ANNAME'][iant]

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Convert uvfits file to filterbank', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='show plots')
    parser.add_argument('--weight', action='store_true', help='Multiply output by weight')
    parser.add_argument('--no-cross', action='store_false', dest='do_cross', default=True)
    parser.add_argument('--maxfiles', help='max files to open', default=np.inf, type=int)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    fin = values.files[0]
    hdu = fits.open(fin)
    g = hdu['PRIMARY']
    hdu.info()
    outfiles = {}
    hdr = g.header
    fcent = hdr['CRVAL4'] / 1e6 # MHz
    foff = hdr['CDELT4'] / 1e6 # MHz
    ref_channel = hdr['CRPIX4'] # pixels

    for row in g.data:
        bl = int(row['BASELINE'])
        ia1 = bl % 256 - 1
        ia2 = bl // 256 - 1
        a1 = get_antname(hdu, ia1)
        a2 = get_antname(hdu, ia2)
        if not (a1 == a2 or values.do_cross):
            continue


        data = row['DATA']
        nchan = data.shape[-3]
        spec = np.zeros(nchan, np.complex64)
        spec.real = data[...,0].reshape(nchan)
        spec.imag = data[...,1].reshape(nchan)
        weights = data[...,2].reshape(nchan)
        if values.weight:
            spec *= weights
            
        if values.show:
            pylab.plot(spec.real)
            pylab.plot(spec.imag)
            pylab.show()

        if bl not in list(outfiles.keys()) and len(outfiles) < values.maxfiles:
            if ia1 == ia2:
                extra = 'auto'
            else:
                extra = 'cross'
                
            outf = '{}-{}-{}-{}.fil'.format(fin.replace('.fits',''),a1,a2, extra)
            jd = row['DATE']
            try:
                inttime = row['INTTIM']
            except:
                inttime = 1
                
            #dayfrac = row['_DATE']
            fulljd = float(jd) 
            mjd = fulljd - 2400000.5
            print('dates', jd,  fulljd, '{:15}'.format(mjd))
            fch1 = fcent - foff*(nchan - ref_channel)
            hdr = {'fch1':fch1, 'foff':foff,'tsamp':inttime, 'tstart':mjd, 'nbits':32, 'nifs':1, 'nchans':nchan, 'src_raj':0.0, 'src_dej':0.0}
            print('Opening', bl, outf, hdr)
            fout = sigproc.SigprocFile(outf, 'w', hdr)
            outfiles[bl] = fout

        if bl in outfiles:
            fout = outfiles[bl]
            abs(spec).tofile(fout.fin)

    for k,fout in outfiles.items():
        fout.fin.close()
        
    
if __name__ == '__main__':
    _main()
