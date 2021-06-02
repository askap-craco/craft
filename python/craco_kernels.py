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
import fdmt
import craco
import boxcar
import uvfits
import craco_plan

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class Kernel(object):
    def __init__(self, uvsource, plan, values):
        self.uvsource = uvsource
        self.plan = plan
        self.values = values

class Prepare(Kernel):
    def __call__(self, dblk):
        dprep = np.zeros((plan.nuvrest, plan.nt, plan.ncin, plan.nuvwide), dtype=np.complex64) # Needs to be zeros otherwise badness?
        
        # prepare data - at this point its just transposing the block
        for irun, fdmtrun in enumerate(plan.fdmt_runs):
            for iuv, uvcell in enmuerate(fdmtrun):
                cell_data = uvd.extract(dblk)
                logging.debug('blkt %s iuv %s cell_data shape=%s', blkt, iuv, cell_data.shape)
                # Truncating times for the moment, as we don't keep history
                dblk[irun, :, :, iuv] = cell_data

        return dprep


class MiniFdmt(Kernel):
    def __call__(self, dblk):
        plan = self.plan

        dfdmt = np.zeros((plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide), dtype=np.complex64)
        for irun, fdmtrun in enumerate(plan.fdmt_plan.fdmt_runs):
            fch1 = min([cell.freqs[0] for cell in fdmtrun]) # effective frequency - should have it as a parameter in plan rather than calculating it
            mincell = min(fdmtrun, key=lambda cell:cell.chan_start)
            assert mincell.freqs[0] == fch1
            minchan = mincell.chan_start
            thefdmt = fdmt.Fdmt(fch1, plan.foff, plan.ncin, plan.ndout, plan.nt) # no history for now
            
            for iuv, uvcell in enumerate(fdmtrun):
                # Truncating times for the moment, as we don't keep history
                tf = dblk.get(uvcell.blid)
                # tf.shape is (nc, nt)
                d = np.zeros((self.plan.ncin, self.plan.nt), dtype=np.complex64)
                subband_start = uvcell.chan_start - minchan
                assert subband_start >= 0
                subband_end = subband_start + uvcell.nchan
                d[subband_start:subband_end, :] = tf[uvcell.chan_slice, :]
                dout = thefdmt(d) # Don't keep history. Shape=(ndout, ndout+nt)
                dfdmt[irun, :, :, iuv] = dout[:, :plan.nt].T
                
        return dfdmt

class Gridder(Kernel):
    '''
    Does Gridding of UV data into a complex grid
    '''
    
    def __call__(self, data1, data2=None):
        '''
        Grids the data
        
        :data1: NUV visibilities
        :data2: If not None, NUV visibilities and the gridder does complex to real gridding,
        so that after FFT you get data1 in the real part and data2 in the imaginary part.
        '''
        npix = self.plan.npix
        g = np.zeros((npix, npix), dtype=self.plan.dtype)
        plan = self.plan.uv_plan
        nuv = len(plan)
        assert data1.shape[0] == nuv, 'UVPlan and grid data have different NUV'
        for iuv in xrange(nuv):
            upix, vpix = map(int, plan[iuv, 2:4])
            v1 = data1[iuv]
            if data2 is not None:
                v2 = data2[iuv]*1j
            else:
                v2 = 0
                
            g[vpix, upix] += v1 + v2
            g[npix-vpix, npix-upix] += np.conj(v1) - np.conj(v2)

        return g

class Imager(Kernel):
    '''
    Takes a grid and makes an image using the FFT
    '''
    def __call__(self, g):
        return image_fft(g)

class Boxcar(Kernel):
    def __init__(self, *args, **kwargs):
        super(Boxcar, self).__init__(*args, **kwargs)
        plan = self.plan
        self.bc = boxcar.ImageBoxcar(plan.nd, plan.npix, plan.nbox, plan.boxcar_weight)
        
    def __call__(self, img):
        return self.bc(img)

class Grouper(Kernel):
    '''
    Groups boxcar candidates
    
    Candidates is a list of (idm, t, xpix, ypix, boxwidth, sn)
    '''
    def __init__(self, *args, **kwargs):
        super(Grouper, self).__init__(*args, **kwargs)        
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

class ImagePipeline(Kernel):
    def __init__(self, *args, **kwargs):
        super(ImagePipeline, self).__init__(*args, **kwargs)
        self.imager = Imager(self.uvsource, self.plan, self.values)
        self.grider = Gridder(self.uvsource, self.plan, self.values)
        self.boxcar = Boxcar(self.uvsource, self.plan, self.values)
        self.grouper = Grouper(self.uvsource, self.plan, self.values)

    def __call__(self):
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
        gridder = self.gridder
        imager = self.imager
        boxcar = self.boxcar
        grouper = self.grouper
        outfname = fname + '.img.dat'
        outgridname= fname + '.grid.dat'
        candname = fname + '.cand'
        npix = self.plan.values.npix
        outshape = (nd, nt/2, npix, npix)
        logging.info("Input shape is %s. Writing output data to %s shape is (nd,nt/2,npix,npix)=%s", d.shape, outfname, outshape)
        fout = open(outfname, 'w')
        gout = open(outgridname, 'w')
        
        assert nt % 2 == 0, 'Nt must be divisible by 2 as were doing complex-to-real gridding'
        
        for idm in xrange(nd):
            for t in xrange(nt/2):
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

class CracoPipeline(Kernel):
    def __init__(self, values):
        uvsource = uvfits.open(values.uv)
        plan = craco_plan.PipelinePlan(uvsource, values)
        super(CracoPipeline, self).__init__(uvsource, plan, values)

        self.image = ImagePipeline(self.uvsource, self.plan, self.values)
        self.fdmt = MiniFdmt(self.uvsource, self.plan, self.values)

    def __call__(self):
        for blkt, d in enumerate(self.uvsource.time_blocks(self.plan.nt)):
            blk = self.fdmt(d)
            dout = self.image(blk)
            return dout


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
