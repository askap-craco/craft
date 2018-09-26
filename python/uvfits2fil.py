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
import sigproc
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
    print g.columns
    outfiles = {}
    hdr = g.header
    fcent = hdr['CRVAL4'] / 1e6 # MHz
    foff = hdr['CDELT4'] / 1e6 # MHz
    ref_channel = hdr['CRPIX4'] # pixels
    fch1 = fcent - foff*ref_channel

    for row in g.data:
        bl = int(row['BASELINE'])
        ia1 = bl % 256 - 1
        ia2 = bl // 256 - 1

        data = row['DATA']
        nchan = data.shape[3]
        spec = np.zeros(nchan, np.complex64)
        spec.real = data[0,0,0,:,0,0]
        spec.imag = data[0,0,0,:,0,1]
        weights = data[0,0,0,:,0,2]
        if values.show:
            pylab.plot(spec.real)
            pylab.plot(spec.imag)
            pylab.show()

        if bl not in outfiles.keys():
            a1 = get_antname(hdu, ia1)
            a2 = get_antname(hdu, ia2)
            if ia1 == ia2:
                extra = 'auto'
            else:
                extra = 'cross'
                
            outf = '{}-{}-{}-{}.fil'.format(fin.replace('.fits',''),a1,a2, extra)
            jd = row['DATE']
            inttime = row['INTTIM']
            dayfrac = row['_DATE']
            fulljd = float(jd) + float(dayfrac)
            mjd = fulljd - 2400000.5
            print 'dates', jd, dayfrac, fulljd, '{:15}'.format(mjd)
            hdr = {'fch1':fch1, 'foff':foff,'tsamp':inttime, 'tstart':mjd, 'nbits':32, 'nifs':1, 'nchans':nchan, 'src_raj':0.0, 'src_dej':0.0}
            print 'Opening', bl, outf, hdr
            fout = sigproc.SigprocFile(outf, 'w', hdr)
            outfiles[bl] = fout

        fout = outfiles[bl]
        abs(spec).tofile(fout.fin)

    for k,fout in outfiles.iteritems():
        fout.fin.close()
        
    
if __name__ == '__main__':
    _main()
