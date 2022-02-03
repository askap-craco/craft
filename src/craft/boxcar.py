#!/usr/bin/env python
"""
Boxcar class - for FDMT pipeline testing

Copyright (C) CSIRO 2020
"""
import numpy as np
import logging

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Boxcar(object):
    def __init__(self, nd, nbox, dtype=np.float32):
        '''
        Make a boxcar object

        :nd: number of DM trials - i.e. number of independant values to boxcar
        :nbox: Number of boxcar trials
        '''
        assert nd > 0, 'Invalid number nd'
        assert nbox > 0, 'Invalid nbox'
        self.history = np.zeros((nd, nbox), dtype=dtype)
        self.nd = nd
        self.nbox = nbox

    def do_boxcar(self, din):
        nd, nt = din.shape
        assert nd == self.nd, 'Invalid input shape={}'.format(din.shape)
        assert nt >= self.nbox, 'Invalid input shape. Need more nt than nobox'
    
        assert din.dtype == self.history.dtype, 'Invalid input dtype'
        dout = np.zeros((nd, nt, self.nbox), dtype=din.dtype)

        # deliberately not efficient
        # We'll just concatenate the history and the input
        # First nbox samples need to use history
        
        dincat = np.hstack([self.history, din])
        nbox = self.nbox
        assert dincat.shape == (self.nd, nt+self.nbox), 'Got funny shape'

        # boxcar 0 == the input
        dout[:, :, 0] = dincat[:, nbox:nt+nbox]

        for b in range(1, nbox):
            dout[:, :, b] = dincat[:, nbox-b:nt+nbox-b] + dout[:, :, b-1]

        self.history[:,:] = din[:, nt-nbox:nt].copy()

        return dout

    __call__ = do_boxcar

class ImageBoxcar(object):
    def __init__(self, nd, npix, nbox, weight='sqrt', dtype=np.float32):
        '''
        Makes an image boxcar - tries to be more faithful to how the the real CRACO image pipeline works

        :nd: Number of DMs to keep a history for
        :npix: Number of pixels on the side of an image
        :nbox: Number of boxcars to compute
        :weight: 'sum' for the sum over the history, 'avg' for the averge over the history or 'sqrt' 
        for the square root over the history. You want 'sqrt' for the S/N weighting.
        '''
        assert nd > 0
        assert npix > 0
        assert nbox > 0
        self.nd = nd
        self.npix = npix
        self.nbox = nbox
        self.history = np.zeros((nd, npix, npix, nbox-1), dtype=dtype)

        if weight == 'sum':
            self.weights = np.ones(nbox, dtype=dtype)
        elif weight == 'avg':
            self.weights = 1./(np.arange(nbox, dtype=dtype) + 1)
        elif weight == 'sqrt':
            self.weights = 1./np.sqrt(np.arange(nbox, dtype=dtype) + 1)
        else:
            raise ValueError('Invalid weight type: %s' % weight)

    def boxcar_image(self, d, img):
        '''
        Does the boxcar on the given image and updates the history
        
        :d: DM index. 0 <= d < nd
        :img: Image shape (npix, npix)
        :returns: (npix, npix, nbox) boxcar output
        '''

        nd = self.nd
        npix = self.npix
        nbox = self.nbox
        assert 0 <= d < nd, 'Invalid dm index'
        assert img.shape == (npix, npix)

        dout = np.zeros((npix, npix, nbox), dtype=self.history.dtype)

        # set initial boxcar value
        dout[:,:, 0] = img

        # Compute boxcar
        for b in range(1, nbox):
            dout[:,:,b] = (dout[:,:,b-1] + self.history[d, :,:,b-1])

        # Multiply by boxcar weights
        dout *= self.weights
            
        # update history - shift everything to the right and put the current image in to position [0]
        self.history[d,:,:, 1:] = self.history[d,:,:,:-1]
        self.history[d,:,:,0] = img

        return dout

    __call__ = boxcar_image

        

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    import pylab

    nd = 10
    npix = 1
    nbox = 6
    
    ib = ImageBoxcar(nd, npix, nbox, 'avg')
    nt = 10

    indata = np.zeros((nt, npix, npix))
    indata[0, :, :] = 1

    outd = []
    for t in range(nt):
        outd.append(ib(0, indata[t,:,:]))

    outd = np.array(outd)

    pylab.plot(outd[:, 0, 0, :])
    pylab.show()
    


    
    

if __name__ == '__main__':
    _main()
