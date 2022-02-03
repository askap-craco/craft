#!/usr/bin/env python
"""
Makes dedispersed visibilities suitable for the CRACO imaging pipeline that does gridding, FFT and boxcar.

New version does partial FDMT with NUVWIDE

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from astropy.io import fits
from scipy import constants
from . import fdmt
import warnings
from .craco import *
from . import craco
from . import uvfits

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def fdmt_baselines(f, uvcells, values):
    freqs = craco.get_freqs(hdul)
    nrows = vis.size
    thefdmt = fdmt.Fdmt(fch1*1e6, foff*1e6, nchan, values.ndm, values.nt)
    nbl = len(baselines)
    baselines = sorted(baselines.keys())
    nuv = len(uvcells)
    for blkt, d in enumerate(time_blocks(vis, values.nt)):
        dblk = np.zeros((nuv, values.ndm, values.nt), dtype=np.complex64)
        alld = []
        for k, v in d.items():
            alld.append(v)
        alld = np.array(alld)
        logging.info('UV data output shape is %s nbl=%s blkt=%s d.shape=%s real/imag mean: %f/%f std: %f/%f',
                     dblk.shape, nbl, blkt, len(d),
                     alld.real.mean(), alld.imag.mean(),
                     alld.real.std(), alld.imag.std())
        # FDMT everything
        for iuv, uvd in enumerate(uvcells):
            cell_data = uvd.extract(d)
            logging.debug('blkt %s iuv %s cell_data shape=%s', blkt, iuv, cell_data.shape)
            # Truncating times for the moment, as we don't keep history
            dblk[iuv, :, :] = thefdmt(cell_data)[:, :values.nt]

            if values.show:
                fig, ax = pylab.subplots(1,2)
                ax[0].imshow(abs(dblk[iuv, :, :]), aspect='auto', origin='lower')
                ax[0].set_xlabel('Time')
                ax[0].set_ylabel('IDM')
                ax[0].set_title('Amplitude iuv={}'.format(iuv))
                ax[1].imshow(np.angle(dblk[iuv, :, :]), aspect='auto', origin='lower')
                ax[1].set_xlabel('Time')
                ax[1].set_ylabel('IDM')
                ax[1].set_title('Phase iuv={}'.format(iuv))
                pylab.show()

        if values.outfile:
            outfile = values.outfile
        else:
            outfile = values.files
            
        fname = '{}.ndm{:d}_nt{:d}.b{:d}.uvdata.{}'.format(outfile, values.ndm, values.nt, blkt, values.format)
        # initial shape is NUV, NDM, NT order
        nuv, ndm, nt = dblk.shape
        dblk *= values.output_scale
        dblk = craco.fdmt_transpose(dblk, ncu=values.nfftcu)

        logging.info('Writing shape %s to %s in format=%s num nonozero=%d scaling with %f real/imag mean:%f/%f std:%f/%f',
                     dblk.shape, fname, values.format, np.sum(dblk != 0), values.output_scale,
                     dblk.real.mean(), dblk.imag.mean(),
                     dblk.real.std(), dblk.imag.std())
        
        if values.format == 'npy':
            np.save(fname, dblk)
        else:
            dblk.tofile(fname)


                    
def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Makes dedispersed visibilities', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--npix', help='Number of pixels in image', type=int, default=256)
    parser.add_argument('--cell', help='Image cell size (arcsec)', default='10,10')
    parser.add_argument('--nt', help='Number of times per block', type=int, default=256)
    parser.add_argument('--ndm', help='Number of DM trials', type=int, default=16)
    parser.add_argument('--outfile', help='Output filename base. Defaults to input filename')
    parser.add_argument('--format', help='Output format', choices=('npy','raw'), default='raw')
    parser.add_argument('--nfftcu', type=int, help='Number of FFT Computing Units for transpose', default=1)
    parser.add_argument('--output-scale', type=float, help='Scalar to multiply output by to change level', default=1.0)
    parser.add_argument('-s','--show', help='Show plots', action='store_true')
                        
    parser.add_argument(dest='files', nargs='?')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Running %s with arguments %s', sys.argv[0], values)
    fdmt_baselines(hdul, baselines, uvcells, values)

if __name__ == '__main__':
    _main()
