#!/usr/bin/env python
"""
Python version of the imaging kernel

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from craco import image_fft

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Gridder(object):
    '''
    Does Gridding of UV data into a complex grid
    '''
    def __init__(self, config, npix, dtype=np.complex64):
        '''
        Creates a gridder.
        :config: UV Grid congirutation = basically a 2D array NUV rows and columns 3,4 containing
        the u,v pixel coordinates 
        :npix: Size of the grid to make
        :dtype: data type of the grid to use
        '''
        self.config = config
        self.npix = npix
        assert 0 < npix
        #assert np.iscomplex(dtype)
        self.dtype = dtype

    
    def __call__(self, data1, data2=None):
        '''
        Grids the data
        
        :data1: NUV visibilities
        :data2: If not None, NUV visibilities and the gridder does complex to real gridding,
        so that after FFT you get data1 in the real part and data2 in the imaginary part.
        '''
        npix = self.npix
        g = np.zeros((npix, npix), dtype=self.dtype)
        nuv = len(self.config)
        assert data1.shape[0] == nuv, 'UVConfig and grid data have different NUV'
        for iuv in xrange(nuv):
            upix, vpix = map(int, self.config[iuv, 2:4])
            v1 = data1[iuv]
            if data2:
                v2 = data2[iuv]*1j
            else:
                v2 = 0
                
            g[vpix, upix] += v1 + v2
            g[-vpix, -upix] += np.conj(v1) - np.conj(v2)

        return g

class Imager(object):
    '''
    Takes a grid and makes an image using the FFT
    '''
    def __call__(self, g):
        return image_fft(g)
        

def image_pipeline(fname, values):
    uvgrid = np.loadtxt(values.uvgrid)
    if fname.endswith('.npy'):
        d = np.load(fname)
        nuv, nd, nt = d.shape
    else:
        nuv = uvgrid.shape[0]
        nd = values.ndm
        nt = values.nt
        d = np.fromfile(fname, dtype=np.complex64).reshape(nuv, nd, nt)

    assert uvgrid.shape[0] == nuv
    assert d.shape == (nuv, nd, nt)
    
    idm = 0
    t = 0
    gridder = Gridder(uvgrid, values.npix)
    imager = Imager()
    outfname = fname + '.img.dat'
    npix = values.npix
    outshape = (nd, nt, npix, npix)
    logging.info("Input shape is %s. Writing outut data to %s shape is %s", d.shape, outfname, outshape)
    fout = open(fname+'.img.raw', 'w')
    
    for idm in xrange(nd):    
        for t in xrange(nt):
            g = gridder(d[:, idm, t])
            img = imager(g)
            img.real.astype(np.complex64).tofile(fout)
            if values.show:
                pylab.imshow(abs(img), aspect='auto', origin='lower')
                pylab.title('idm={} t={}'.format(idm, t))
                pylab.show()

    logging.info("Wrote outpuut images to %s shape=", outfname, outshape)

    fout.close()
            

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--uvgrid', help='UV Grid file output')
    parser.add_argument('--npix', help='Number of pixels', type=int, default=256)
    parser.add_argument('--outfile', help='Raw data output file')
    parser.add_argument('--nt', type=int, help='Number of times per block. Used if using raw format', default=256)
    parser.add_argument('--ndm', type=int, help='Number of DMs.  Used if raw format', default=16)
    parser.add_argument('-s','--show', action='store_true', help='Show plots', default=False)
    parser.add_argument(dest='files', nargs='*')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for f in values.files:
        image_pipeline(f, values)
    
    

if __name__ == '__main__':
    _main()
