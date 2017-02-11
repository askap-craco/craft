#!/usr/bin/env python
"""
Plot a psrdada file

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'

import logging
import dada
import pylab
import numpy as np
from cmdline import strrange

class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.01f}   '.format(x, y, z)


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-t', '--time', help='Sample times', default='0,256')
    parser.add_argument('-b', '--bmax', help='Maximum beam to plot', type=int, default=2)
    parser.add_argument('--imrange', help='Imavge vertical plot range')
    parser.add_argument('--nxy', help='nxy', default='6,12')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    f = dada.DadaFile(values.files[0])
    nbeam = int(f.hdr['NBEAM'])
    nchan = int(f.hdr['NCHAN'])
    npol = int(f.hdr['NPOL'])
    tstart, nint = map(int, values.time.split(','))
    dtype = np.dtype(f.hdr.get('DTYPE', '<f4'))
    order = f.hdr.get('DORDER', 'TFBP')
    sz = dtype.itemsize
    byteoff = 4*nchan*nbeam*npol*tstart
    f.fin.seek(byteoff + f.header_size)
    while True:
        plot(f, dtype, nint, nbeam, nchan, npol, order, values)


def plot(f, dtype, nint, nbeam, nchan, npol, order, values):
    nelements = nbeam*nchan*nint*npol

    v = np.fromfile(f.fin, dtype=dtype, count=nelements)
    
    if order == 'TBPF':
        v.shape = (nint, nbeam, npol, nchan)
    elif order == 'TFBP':
        v.shape = (nint, nchan, nbeam, npol)
    else:
        raise ValueError('Unknown order %s' % order)


    print 'Got dtypes', dtype, 'ORDER', order, 'shape', v.shape

    fig, axes = pylab.subplots(*map(int, values.nxy.split(',')))
    if values.imrange:
        vmin, vmax = map(float, values.imrange.split(','))
    else:
        vmin = None
        vmax = None
    
    for iax, ax in enumerate(axes.flat):
        assert order=='TFBP'
        pol = iax / nbeam
        beam = iax % nbeam
        img = v[:, :, beam,pol].T
        print 'beam', beam, 'pol', pol,'max/min/mean/rms {}/{}/{}/{}'.format(img.max(), img.min(), img.mean(), img.std())
        im = ax.imshow(img, aspect='auto', vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
        ax.text(0, 0, 'beam %d pol %d' % (beam, pol), va='top', ha='left') 
        ax.format_coord = Formatter(im)


    pylab.subplots_adjust(wspace=0, hspace=0)
    #pylab.tight_layout()
    pylab.show()
    

if __name__ == '__main__':
    _main()
