#!/usr/bin/env python
"""
Plot a psrdada file

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'

import logging
from . import dada
import pylab
import numpy as np
from .cmdline import strrange
from .crafthdr import DadaHeader

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
    parser.add_argument('--tscrunch', help='Tscrunch factor', type=int, default=1)
    parser.add_argument('--fscrunch', help='Fscrunch factor', type=int, default=1)
    parser.add_argument('--imrange', help='Imavge vertical plot range')
    parser.add_argument('--nxy', help='nxy', default='6,12')
    parser.add_argument('--order', help='Force ordering')
    parser.add_argument('--mask-limit', help='Mask values with absolute less than this', default=0.0, type=float)
    parser.add_argument('-r','--rescale', action='store_true', help='Rescale data', default=False)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    #f = DadaFile(values.files[0])
    hdr = DadaHeader.fromfile(values.files[0])
    nbeam = int(hdr['NBEAM'][0])
    nchan = int(hdr['NCHAN'][0])
    npol = int(hdr['NPOL'][0])
    tstart, nint = list(map(int, values.time.split(',')))
    dtype = np.dtype(hdr.get_value('DTYPE', '<f4'))
    if values.order is not None:
        order = values.order
    else:
        order = hdr.get_value('DORDER', 'TFBP')
        
    nint = int(hdr.get_value('NT', nint))
    transpose = None



    if order == 'TFBP':
        shape = (nint, nchan, nbeam, npol)
    elif order == 'TBPF':
        shape = (nint, nbeam, npol, nchan)
        transpose = [0, 3, 1, 2]
    elif order == 'BPFT':
        shape = (nbeam, npol, nchan, nint)
        transpose = [3, 2, 0, 1]
    else:
        raise ValueError('Unknown order {}'.format(order))

    print(nbeam, nchan, npol, tstart, nint, dtype, order, nint, shape, transpose)
    
    f = dada.DadaFile(values.files[0], shape=shape)
    for b in f.blocks():
        orig_shape = b.shape
        if transpose is not None:
            b = b.transpose(transpose)

        b = np.ma.masked_where(abs(b) <= values.mask_limit, b)

        # At this point b.shape is (Time, Chan, Beam, Pol)
        print(b.shape, shape, transpose, orig_shape)
        tscrunch = values.tscrunch
        fscrunch = values.fscrunch
        nint_out = nint/tscrunch
        nchan_out = nchan/fscrunch
        b = b.reshape(nint_out, tscrunch, nchan, nbeam, npol).mean(axis=1)
        b = b.reshape(nint_out, nchan_out, fscrunch, nbeam, npol).mean(axis=2)
        print(b.shape, shape, transpose, orig_shape)
        plot(b, nint_out, nbeam, nchan_out, npol, values)


def plot(v, nint, nbeam, nchan, npol, values):
    fig, axes = pylab.subplots(*list(map(int, values.nxy.split(','))), sharex=True, sharey=True)
    if values.imrange:
        vmin, vmax = list(map(float, values.imrange.split(',')))
    else:
        vmin = None
        vmax = None

    nbeampol = nbeam*npol
    for iax, ax in enumerate(axes.flatten()[:nbeampol]):
        assert v.shape == (nint, nchan, nbeam, npol), 'Invalid shape {} expected ({},{},{},{})'.format(v.shape, nint, nchan, nbeam, npol) # assume ordering TFBP
        pol = iax / nbeam
        beam = iax % nbeam
        if pol > npol or beam > nbeam:
            break
        
        img = v[:, :, beam,pol]

        if values.rescale:
            img -= img.mean(axis=0)
            img /= img.std(axis=0)
            print('RESCALED beam', beam, 'pol', pol,'max/min/mean/rms {}/{}/{}/{}'.format(img.max(), img.min(), img.mean(), img.std()), img.shape)
            
        im = ax.imshow(img.T, aspect='auto', vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
        ax.text(0, 0, 'beam %d pol %d' % (beam, pol), va='top', ha='left') 
        ax.format_coord = Formatter(im)


    pylab.subplots_adjust(wspace=0, hspace=0)
    #pylab.tight_layout()
    pylab.show()
    

if __name__ == '__main__':
    _main()
