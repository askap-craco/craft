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
from . import craco
from .craco import image_fft, printstats
from .craco_kernels import Imager
from .boxcar import ImageBoxcar

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
        for iuv in range(nuv):
            upix, vpix = list(map(int, self.config[iuv, 2:4]))
            v1 = data1[iuv]
            if data2 is not None:
                v2 = data2[iuv]*1j
            else:
                v2 = 0
                
            g[vpix, upix] += v1 + v2
            g[npix-vpix, npix-upix] += np.conj(v1) - np.conj(v2)

        return g

class Grouper(object):
    '''
    Groups boxcar candidates
    
    Candidates is a list of (idm, t, xpix, ypix, boxwidth, sn)
    '''
    def __init__(self, threshold):
        self.threshold = threshold
        self.candidates = []

    def __call__(self, idm, t, boxout):
        '''
        runs on a single boxcar output with shape (npix, npix, nbox) 
        and returns the associated candidates
        :idm: DM value
        :t: Time stamp
        :boxout: boxcar output shape (npix, npix, nbox)
        '''

        # Find the boxcar with the largest value for each pixel
        # See https://stackoverflow.com/questions/42519475/python-numpy-argmax-to-max-in-multidimensional-array
        max_box = np.argmax(boxout, axis=2)
        i,j = np.indices(max_box.shape)
        
        # Find pixels where the largest boxcar is larger than the threshold
        pys, pxs = np.where(boxout[i,j,max_box] >= self.threshold)

        cands = []
        for xpix, ypix in zip(pxs, pys):
            boxwidth = max_box[ypix, xpix] # boxwidth = 1 for the smallest boxcar
            sn = boxout[ypix,xpix,boxwidth]
            cand = (idm, t, xpix, ypix, boxwidth+1, sn)
            cands.append(cand)

        self.candidates.extend(cands)
        
        return cands

    def to_file(self, fname):
        logging.info('Saved %s candiates to %s', len(self.candidates), fname)
        fmt = '%d %d %d %d %d %0.1f'
        np.savetxt(fname, np.array(self.candidates), header='idm t xpix ypix boxwidth sn', fmt=fmt)

def image_pipeline(fname, values):
    uvgrid = np.loadtxt(values.uvgrid)
    if uvgrid.ndim == 1:
        uvgrid = uvgrid[np.newaxis, :]
        
    if fname.endswith('.npy'):
        d = np.load(fname)
        ncu, nd, nt_on_ncu, nuv = d.shape
        nt = nt_on_ncu * ncu
        # Output expected to be (nd, nt, nuv)
        d = np.transpose(d, (1, 2, 0, 3)).reshape(nd, nt, nuv)
    else:
        nuv = uvgrid.shape[0]
        nd = values.ndm
        nt = values.nt
        ncu = values.nfftcu
        d = np.fromfile(fname, dtype=np.complex64)
        d = craco.fdmt_transpose_inv(d, ncu=values.nfftcu, ndm=nd, nt=nt, nuv=nuv)
        # Image transpose outputs (nuv, ndm, nt) - should fix everything to be consistent
        # But for now we'll transpose to (nd, nt, nuv) as this is what the next code expects
        d = np.transpose(d, (1, 2, 0))

    assert uvgrid.shape[0] == nuv
    assert d.shape == (nd, nt, nuv)
    
    idm = 0
    t = 0
    gridder = Gridder(uvgrid, values.npix)
    imager = Imager()
    boxcar = ImageBoxcar(nd, values.npix, values.nbox, values.boxcar_weight)
    grouper = Grouper(values.threshold)
    outfname = fname + '.img.dat'
    outgridname= fname + '.grid.dat'
    candname = fname + '.cand'
    npix = values.npix
    outshape = (nd, nt/2, npix, npix)
    logging.info("Input shape is %s. Writing output data to %s shape is (nd,nt/2,npix,npix)=%s", d.shape, outfname, outshape)
    fout = open(outfname, 'w')
    gout = open(outgridname, 'w')

    assert nt % 2 == 0, 'Nt must be divisible by 2 as were doing complex-to-real gridding'

    for idm in range(nd):
        for t in range(nt/2):
            g = gridder(d[idm, 2*t, :], d[idm, 2*t+1, :])
            g.tofile(gout)
            img = imager(g).astype(np.complex64)
            img.tofile(fout)
            grouper(idm, 2*t, boxcar(idm, img.real))
            grouper(idm, 2*t + 1, boxcar(idm, img.imag))

            rlabel = 'img real idm={} t={}'.format(idm, t)
            ilabel = 'img imag idm={} t={}'.format(idm, t+1)
            printstats(img.real, rlabel)
            printstats(img.imag, ilabel)
            printstats(g.real, 'grid.real')
            printstats(g.imag, 'grid.imag')
            if values.show:
                fig, ax = pylab.subplots(2,2)
                ax[0,0].imshow(img.real, aspect='auto', origin='lower')
                ax[0,1].imshow(img.imag, aspect='auto', origin='lower')
                ax[1,0].imshow(g.real, aspect='auto', origin='lower')
                ax[1,1].imshow(g.imag, aspect='auto', origin='lower')
                ax[0,0].set_title(rlabel)
                ax[0,1].set_title(ilabel)
                ax[1,0].set_title('real(UV plane)')
                ax[1,1].set_title('imag(UV plane)')
                pylab.show()

    logging.info("Wrote output images to %s shape=%s (nd,nt,npix,npix)=dtype=%s", outfname, outshape, img.dtype)

    grouper.to_file(candname)
    fout.close()
    gout.close()
            

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--uvgrid', help='UV Grid file output')
    parser.add_argument('--npix', help='Number of pixels', type=int, default=256)
    parser.add_argument('--outfile', help='Raw data output file')
    parser.add_argument('--nt', type=int, help='Number of times per block. Used if using raw format', default=256)
    parser.add_argument('--ndm', type=int, help='Number of DMs.  Used if raw format', default=16)
    parser.add_argument('--nfftcu', type=int, help='Number of FFT Computing Units for transpose', default=1)
    parser.add_argument('--nbox', type=int, help='Number of boxcars to compute', default=8)
    parser.add_argument('--threshold', type=float, help='Threshold for candidate grouper', default=10)
    parser.add_argument('--boxcar-weight', choices=('sqrt','avg','sum'), help='Boxcar weight type', default='sqrt')
    parser.add_argument('-s','--show', action='store_true', help='Show plots', default=False)
    parser.add_argument(dest='files', nargs='*')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Running %s with arguments %s', sys.argv[0], values)
    for f in values.files:
        image_pipeline(f, values)
    
    

if __name__ == '__main__':
    _main()
