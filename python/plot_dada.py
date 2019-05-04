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
from crafthdr import DadaHeader

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
    tstart, nint = map(int, values.time.split(','))
    dtype = np.dtype(hdr.get_value('DTYPE', '<f4'))
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


    f = dada.DadaFile(values.files[0], shape=shape)
    for b in f.blocks():
        plot(b, nint, nbeam, nchan, npol, transpose, values)


def plot(v, nint, nbeam, nchan, npol, transpose, values):

    if transpose is not None:
        v.transpose(transpose)


    fig, axes = pylab.subplots(*map(int, values.nxy.split(',')), sharex=True, sharey=True)
    if values.imrange:
        vmin, vmax = map(float, values.imrange.split(','))
    else:
        vmin = None
        vmax = None
    
    for iax, ax in enumerate(axes.flat):
        assert v.shape == (nint, nchan, nbeam, npol) # assume ordering TFBP
        pol = iax / nbeam
        beam = iax % nbeam
        if pol > npol or beam > nbeam:
            break
        
        img = v[:, :, beam,pol]
        print 'beam', beam, 'pol', pol,'max/min/mean/rms {}/{}/{}/{}'.format(img.max(), img.min(), img.mean(), img.std()), img.shape
        if values.rescale:
            img -= img.mean(axis=0)
            img /= img.std(axis=0)
            print 'RESCALED beam', beam, 'pol', pol,'max/min/mean/rms {}/{}/{}/{}'.format(img.max(), img.min(), img.mean(), img.std()), img.shape
            
        im = ax.imshow(img.T, aspect='auto', vmin=vmin, vmax=vmax, interpolation='none', origin='lower')
        ax.text(0, 0, 'beam %d pol %d' % (beam, pol), va='top', ha='left') 
        ax.format_coord = Formatter(im)


    pylab.subplots_adjust(wspace=0, hspace=0)
    #pylab.tight_layout()
    pylab.show()
    

if __name__ == '__main__':
    _main()
