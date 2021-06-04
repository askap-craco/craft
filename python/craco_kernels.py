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
from craco import printstats

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def imshow_complex(axs, d, title=''):
    axs[0].imshow(d.real, aspect='auto', origin='lower')
    axs[0].set_title(title + ' real')
    axs[1].imshow(d.imag, aspect='auto', origin='lower')
    axs[1].set_title(title + ' imag')

    #logging.debug('%s real %s', title, printstats(d.real))
    #logging.debug('%s imag %s', title, printstats(d.imag))

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
                #logging.debug('irun=%d iuv=%d blid=%s real=%s', irun, iuv, uvcell.blid, craco.printstats(tf.real))
                # tf.shape is (nc, nt)
                d = np.zeros((self.plan.ncin, self.plan.nt), dtype=np.complex64)
                subband_start = uvcell.chan_start - minchan
                assert subband_start >= 0
                subband_end = subband_start + uvcell.nchan
                d[subband_start:subband_end, :] = tf[uvcell.chan_slice, :]
                dout = thefdmt(d) # Don't keep history. Shape=(ndout, ndout+nt)

                if self.plan.values.show_fdmt:
                    fig, axs = pylab.subplots(3,2)
                    imshow_complex(axs[0,:], tf, 'FDMT TF')
                    imshow_complex(axs[1,:], d, 'FDMT d')
                    imshow_complex(axs[2,:], dout, 'FDMT dout')
                    
                    #fig.set_title('irun={} iuv={} blid={}'.format(irun, iuv, uvcell.blid))
                    pylab.show()
                    
                dfdmt[irun, :, :, iuv] = dout[:, :plan.nt].T

        # rescale for ... goodness
        scale_factor = self.plan.fdmt_scale / float(self.plan.nf) # the 2 I think happens becuase of the complex->real gridding

        dfdmt *= scale_factor
        
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

def idm_cff(fch1, plan):
    '''
    Calculate the CFF coefficient to calculate the requred IDM inside an FDMT
    Multiply by the desired IDT and round to get an index


    :idt: Overall IDT we're trying to calculate
    :fch1: Frequency of the bottom of the FDMT channels
    :plan: craco plan
    '''
    fdmt_band = plan.ncin*plan.foff # Bandwidth of FDMT
    f1 = fch1
    f2 = f1 + fdmt_band
    fmin = plan.fmin
    fmax = plan.fmax
    dmcff = fdmt.cff(f2, f1, fmax, fmin)
    assert dmcff >= 0
    return dmcff

def offset_cff(fch1, plan):
    '''
    Calculate the CFF coefficient to calculate the required OFFSET inside an FDMT
    Multiply by the desired IDT and round to get an index

    :idt: Overall IDT we're trying to calculate
    :fch1: Frequency of the bottom of the FDMT channels
    :plan: craco plan
    '''
    f1 = fch1
    fmin = plan.fmin
    fmax = plan.fmax
    ocff = fdmt.cff(f1, fmin, fmax, fmin) + 0.5/float(plan.nd) # Need to add a bit extra. not sure why.
    assert ocff >= 0
    return ocff

class FdmtGridder(Kernel):
    '''
    Grids the data assuming the FDMT has only done some of the dedispersion
    and we need do do the rest of the processing
    blk shape should be  dfdmt = np.zeros((plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide), dtype=np.complex64)
    Also sums into uv pixels which are duplicated
    Returns a uvgrid size (npix, npix) dtype=complex64
    '''
    
    def __call__(self, idm, t, blk):
        assert idm < self.plan.nd
        assert 2*t+1 < self.plan.nt
        plan = self.plan
        #g = gridder(d[idm, 2*t, :], d[idm, 2*t+1, :])
        assert blk.shape == (plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide)
        npix = self.plan.npix
        g = np.zeros((npix, npix), dtype=self.plan.dtype)
        fdmt_band = plan.ncin*plan.foff # Bandwidth of FDMT

        for irun, fdmtrun in enumerate(plan.fdmt_plan.fdmt_runs):
            mincell = min(fdmtrun, key=lambda cell:cell.chan_start)
            minchan = mincell.chan_start
            fch1 = mincell.freqs[0]
            blkdm = int(np.round(idm*idm_cff(fch1, plan)))
            toff = int(np.round(idm*offset_cff(fch1, plan)))
            blkt = 2*t - toff
            logging.debug('Fdmtgrider idm=%d t=%d irun=%d fch1=%0.1f blkdm=%d toff=%d blkt=%d', idm, t, irun, fch1, blkdm, toff, blkt)
            for iuv, uvcell in enumerate(fdmtrun):
                upix, vpix = uvcell.uvpix
                if (blkt >= 0):
                    v1 = blk[irun, blkt+0, blkdm, iuv]
                else:
                    v1 = complex(0,0)

                if (blkt+1 >= 0):
                    v2 = blk[irun, blkt+1, blkdm, iuv]
                else:
                    v2 = complex(0,0)

                # This is summing because some UV Cells get split to do the FDMT and we need to recombine them
                g[vpix, upix] += v1 + v2
                g[npix-vpix, npix-upix] += np.conj(v1) - np.conj(v2)
                
        return g


class Imager(Kernel):
    '''
    Takes a grid and makes an image using the FFT
    '''
    def __call__(self, g):
        img = craco.image_fft(g)
        scale_factor = self.plan.fft_scale/(self.plan.nbl * 2.0) # factor of 2 is becauseo complex-to-real.
        img *= scale_factor
        return img

class Boxcar(Kernel):
    def __init__(self, *args, **kwargs):
        super(Boxcar, self).__init__(*args, **kwargs)
        plan = self.plan
        self.bc = boxcar.ImageBoxcar(plan.nd, plan.npix, plan.nbox, plan.boxcar_weight)
        
    def __call__(self, idm, img):
        return self.bc(idm, img)

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
        threshold = self.plan.threshold
        
        # Find pixels where the largest boxcar is larger than the threshold
        pys, pxs = np.where(boxout[i,j,max_box] >= threshold)

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
        self.gridder = FdmtGridder(self.uvsource, self.plan, self.values)
        self.boxcar = Boxcar(self.uvsource, self.plan, self.values)
        self.grouper = Grouper(self.uvsource, self.plan, self.values)

    def load_from_file(self, fname):
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

        return d

    def __call__(self, blk):
        # blk shape should be  dfdmt = np.zeros((plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide), dtype=np.complex64)
        plan = self.plan
        assert blk.shape == (plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide)
        gridder = self.gridder
        imager = self.imager
        boxcar = self.boxcar
        grouper = self.grouper
        fname = 'testing'
        outfname = fname + '.img.dat'
        outgridname= fname + '.grid.dat'
        candname = fname + '.cand'
        npix = self.plan.npix
        nd = self.plan.nd
        nt = self.plan.nt
        outshape = (nd, nt/2, npix, npix)
        logging.info("Input shape is %s. Writing output data to %s shape is (nd,nt/2,npix,npix)=%s", blk.shape, outfname, outshape)
        fout = open(outfname, 'w')
        gout = open(outgridname, 'w')
        
        assert nt % 2 == 0, 'Nt must be divisible by 2 as were doing complex-to-real gridding'
        
        for idm in xrange(nd):
            for t in xrange(nt/2):
                #g = gridder(d[idm, 2*t, :], d[idm, 2*t+1, :])
                g = gridder(idm, t, blk)
                #g.tofile(gout)
                img = imager(g).astype(np.complex64)
                #img.tofile(fout)
                grouper(idm, 2*t, boxcar(idm, img.real))
                grouper(idm, 2*t + 1, boxcar(idm, img.imag))
                rlabel = 'imag idm={} t={}'.format(idm, t)
                ilabel = 'real idm={} t={}'.format(idm, t)
                logging.debug('img.real idm=%d t=%d %s', idm, t, printstats(img.real, rlabel))
                logging.debug('img.imag idm=%d t=%d %s', idm, t, printstats(img.imag, ilabel))
                logging.debug('grid.real idm=%d t=%d %s', idm, t, printstats(g.real, 'grid.real'))
                logging.debug('grid.imag idm=%d t=%d %s', idm, t, printstats(g.imag, 'grid.imag'))
                if self.values.show_image:
                    fig, ax = pylab.subplots(2,2)
                    imshow_complex(ax[0,:], img, 'image idm={} t={}'.format(idm, t))
                    imshow_complex(ax[1,:], g, 'grid idm={} t={}'.format(idm, t))
                    pylab.show()

            logging.info("Wrote output images to %s shape=%s (nd,nt,npix,npix)=dtype=%s", outfname, outshape, img.dtype)
    
        grouper.to_file(candname)
        fout.close()
        gout.close()

def savefile(fname, arr):
    logging.info('Saving file %s shape=%s dtype=%s', fname, arr.shape, arr.dtype)
    np.save(fname, arr)

class CracoPipeline(Kernel):
    def __init__(self, values):
        uvsource = uvfits.open(values.uv)
        plan = craco_plan.PipelinePlan(uvsource, values)
        super(CracoPipeline, self).__init__(uvsource, plan, values)
        self.fdmt = MiniFdmt(self.uvsource, self.plan, self.values)
        self.image = ImagePipeline(self.uvsource, self.plan, self.values)

    def __call__(self):
        for blkt, d in enumerate(self.uvsource.time_blocks(self.plan.nt)):
            if self.plan.values.save:
                savefile('blk_{}_input.npy'.format(blkt), craco.bl2array(d))

            blk = self.fdmt(d)

            if self.plan.values.save:
                savefile('blk_{}_fdmt.npy'.format(blkt), blk)

            dout = self.image(blk)

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
