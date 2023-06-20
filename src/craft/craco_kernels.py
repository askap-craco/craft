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

from . import fdmt
from . import craco
from . import boxcar
from . import uvfits
from . import craco_plan

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def imshow_complex(axs, d, title=''):
    axs[0].imshow(d.real, aspect='auto', origin='lower')
    axs[0].set_title(title + ' real')
    axs[1].imshow(d.imag, aspect='auto', origin='lower')
    axs[1].set_title(title + ' imag')

    #logging.debug('%s real %s', title, craco.printstats(d.real))
    #logging.debug('%s imag %s', title, craco.printstats(d.imag))

class Kernel(object):
    def __init__(self, uvsource, plan, values):
        self.uvsource = uvsource
        self.plan = plan
        self.values = values

class Prepare(Kernel):
    def __call__(self, dblk):
        plan = self.plan
        dprep = np.zeros((plan.nuvrest, plan.nt, plan.ncin, plan.nuvwide), dtype=np.complex64) # Needs to be zeros otherwise badness?
        
        # prepare data - at this point its just transposing the block
        for irun, fdmtrun in enumerate(plan.fdmt_plan.fdmt_runs):
            minchan = plan.fdmt_plan.run_chan_starts[irun]
            for iuv, uvcell in enumerate(fdmtrun):
                tf = dblk.get(uvcell.blid) # lookup in dictionary - shape=(nchan, nt)
                subband_start = uvcell.chan_start - minchan
                assert subband_start >= 0
                subband_end = subband_start + uvcell.nchan
                assert subband_end <= self.plan.ncin
                subband_slice = slice(subband_start, subband_end)
                dprep[irun, :, subband_slice, iuv] = tf[uvcell.chan_slice, :].T

        return dprep


class MiniFdmt(Kernel):
    def __call__(self, dblk):
        plan = self.plan
        assert dblk.shape == (plan.nuvrest, plan.nt, plan.ncin, plan.nuvwide)

        dfdmt = np.zeros((plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide), dtype=np.complex64)
        for irun, fdmtrun in enumerate(plan.fdmt_plan.fdmt_runs):
            fch1 = self.plan.fdmt_plan.run_fch1[irun] # lookup lookup table
            thefdmt = fdmt.Fdmt(fch1, plan.foff, plan.ncin, plan.ndout, plan.nt) # no history for now
            
            for iuv, uvcell in enumerate(fdmtrun):
                # Truncating times for the moment, as we don't keep history
                d = dblk[irun, :, :, iuv]
                dout = thefdmt(d.T).T # Don't keep history. Shape=(ndout, ndout+nt)

                if self.plan.values.show_fdmt:
                    fig, axs = pylab.subplots(2,2)
                    imshow_complex(axs[0,:], d.T, 'FDMT d')
                    imshow_complex(axs[1,:], dout.T, 'FDMT dout')
                    
                    #fig.set_title('irun={} iuv={} blid={}'.format(irun, iuv, uvcell.blid))
                    pylab.show()
                    
                dfdmt[irun, :, :, iuv] = dout[:plan.nt, :]

        # rescale for ... goodness
        scale_factor = self.plan.fdmt_scale / float(self.plan.nf) # the 2 I think happens becuase of the complex->real gridding

        dfdmt *= scale_factor
        
        return dfdmt

class Gridder(Kernel):
    '''
    Does Gridding of UV data into a complex grid
    '''
    
    def __call__(self, block):
        '''
        Performs complex to complex gridding of the visibility data
        
        Params
        ------
        block: np.ndarray or numpy.ma.core.MaskedArray or dict
                Block containing [nbl, nf, nt] complex visibility data (if array)
                Dict containing [nf, nt] complex visibility data for nbl baselines (if dict)
        '''
        if type(block) not in [np.ndarray, np.ma.core.MaskedArray, dict]:
            raise Exception(f"I expected a np.ndarray/masked_array/dict, but got {type(block)}")
        
        if type(block) == dict:
            expected_nbl = self.plan.nbl
            assert len(block) == expected_nbl, "Vis dict did not have expected no. of baselines"
            blids = list(block.keys())
            data0 = block[blids[0]]
            assert data0.ndim == 2, f"Expected a 2-D array for each baseline, got-{data0.ndim}"
            assert data0.shape[0] == self.plan.nf
            nt = data0.shape[-1]

        else:
            assert block.ndim == 3, "block needs to have [nbl, nf, nt] shape"
            assert block.shape[:2] == (self.plan.nbl, self.plan.nf)
            nt = block.shape[-1]

        npix = self.plan.npix
        assert nt >= 2, "Block needs to have atleast 2 time samples to perform CPLX to CPLX gridding"
        assert nt%2 == 0, "Block needs to have even number of samples"
        g = np.zeros((npix, npix, nt//2), dtype=self.plan.dtype)
        nuv = len(self.plan.uvcells)

        for iuv in range(nuv):
            cell = self.plan.uvcells[iuv]
            upix, vpix = cell.uvpix

            if type(block) == dict:
                v1 = block[cell.blid][cell.chan_slice, ::2].sum(axis=0)
                v2 = block[cell.blid][cell.chan_slice, 1::2].sum(axis=0) * 1j
            else:
                bl_idx = np.where(self.plan.baseline_order == cell.blid)[0][0]
                v1 = block[bl_idx, cell.chan_slice, ::2].sum(axis=0)
                v2 = block[bl_idx, cell.chan_slice, 1::2].sum(axis=0) * 1j
                
            g[vpix, upix, :] += v1 + v2
            g[-vpix, -upix, :] += np.conj(v1) - np.conj(v2)

        return g
    
    def complx_to_real_gridder(self, block):
        '''
        Performs complex to real gridding of the visibility data
        
        Params
        ------
        block: np.ndarray or numpy.ma.core.MaskedArray or dict
                Block containing [nbl, nf, nt] complex visibility data (if array)
                Dict containing [nf, nt] complex visibility data for nbl baselines (if dict)
        '''
        if type(block) not in [np.ndarray, np.ma.core.MaskedArray, dict]:
            raise Exception(f"I expected a np.ndarray/masked_array/dict, but got {type(block)}")
        
        if type(block) == dict:
            expected_nbl = self.plan.nbl
            assert len(block) == expected_nbl, "Vis dict did not have expected no. of baselines"
            blids = list(block.keys())
            data0 = block[blids[0]]
            assert data0.ndim == 2, f"Expected a 2-D array for each baseline, got-{data0.ndim}"
            assert data0.shape[0] == self.plan.nf
            nt = data0.shape[-1]

        else:
            assert block.ndim == 3, "block needs to have [nbl, nf, nt] shape"
            assert block.shape[:2] == (self.plan.nbl, self.plan.nf)
            nt = block.shape[-1]

        npix = self.plan.npix
        assert nt >= 2, "Block needs to have atleast 2 time samples to perform CPLX to CPLX gridding"
        assert nt%2 == 0, "Block needs to have even number of samples"
        g = np.zeros((npix, npix, nt), dtype=self.plan.dtype)
        nuv = len(self.plan.uvcells)

        for iuv in range(nuv):
            cell = self.plan.uvcells[iuv]
            upix, vpix = cell.uvpix

            if type(block) == dict:
                v1 = block[cell.blid][cell.chan_slice,:].sum(axis=0)
                v2 = 0
            else:
                bl_idx = np.where(self.plan.baseline_order == cell.blid)[0][0]
                v1 = block[bl_idx, cell.chan_slice, :].sum(axis=0)
                v2 = 0
                
            g[vpix, upix, :] += v1 + v2
            g[-vpix, -upix, :] += np.conj(v1) - np.conj(v2)

        return g
    
    def grid_with_uvws(self, block, uvws):
        '''
        Performs complex to complex gridding of the visibility data, taking into account the changing uvw of each time sample
        
        Params
        ------
        block: np.ndarray or numpy.ma.core.MaskedArray or dict
                Block containing [nbl, nf, nt] complex visibility data (if array)
                Dict containing [nf, nt] complex visibility data (keyed by blid) for nbl baselines (if dict)
        uvws: dict
                Dict containing UVW values (keyed by blid) as an numpy array of shape [3, nt]
        '''
        if type(block) not in [np.ndarray, np.ma.core.MaskedArray, dict]:
            raise Exception(f"I expected a np.ndarray/masked_array/dict, but got {type(block)}")
        
        if type(block) == dict:
            expected_nbl = self.plan.nbl
            assert len(block) == expected_nbl, "Vis dict did not have expected no. of baselines"
            blids = list(block.keys())
            data0 = block[blids[0]]
            assert data0.ndim == 2, f"Expected a 2-D array for each baseline, got-{data0.ndim}"
            assert data0.shape[0] == self.plan.nf
            nt = data0.shape[-1]

        else:
            assert block.ndim == 3, "block needs to have [nbl, nf, nt] shape"
            assert block.shape[:2] == (self.plan.nbl, self.plan.nf)
            nt = block.shape[-1]

        if type(uvws) != dict:
            raise Exception(f"I expected UVWs to be passed as a dict, but got {type(uvws)}")
        
        first_blid = list(uvws.keys())[0]
        uvws_shape = uvws[first_blid].shape
        assert uvws_shape == (3, nt), f"UVWs shape needs to be (3, {nt}), got {uvws_shape}"

        npix = self.plan.npix
        assert nt >= 2, "Block needs to have atleast 2 time samples to perform CPLX to CPLX gridding"
        assert nt%2 == 0, "Block needs to have even number of samples"
        
        g = np.zeros((npix, npix, nt//2), dtype=self.plan.dtype)
        for t in range(nt//2):
            this_uvw = {}
            for blid, uvw_data in list(uvws.items()):
                this_uvw[blid] = np.array(tuple(uvw_data[:, t]), dtype=[('UU', 'f8'), ('VV', 'f8'), ('WW', 'f8')])
            current_uvcells = craco_plan.get_uvcells(baselines=this_uvw, uvcell=self.plan.uvcell, freqs = self.plan.freqs, Npix = self.values.npix)
            nuv = len(current_uvcells)

            for iuv in range(nuv):
                cell = current_uvcells[iuv]
                upix, vpix = cell.uvpix

                if type(block) == dict:
                    v1 = block[cell.blid][cell.chan_slice, t].sum(axis=0)
                    v2 = block[cell.blid][cell.chan_slice, t+1].sum(axis=0) * 1j
                else:
                    bl_idx = np.where(self.plan.baseline_order == cell.blid)[0][0]
                    v1 = block[bl_idx, cell.chan_slice, t].sum(axis=0)
                    v2 = block[bl_idx, cell.chan_slice, t+1].sum(axis=0) * 1j

                g[vpix, upix, t] += v1 + v2
                g[-vpix, -upix, t] += np.conj(v1) - np.conj(v2)

        return g

def idm_cff(fch1, plan):
    '''
    Calculate the CFF coefficient to calculate the requred IDM inside an FDMT
    Multiply by the desired IDT and round to get an index


    :idt: Overall IDT we're trying to calculate
    :fch1: Frequency of the bottom of the FDMT channels
    :plan: craco plan
    '''
    fdmt_band = plan.ncin*plan.foff # Bandwidth of sub-FDMT
    fmin = plan.fmin # center freq of bottom channel
    fmax = plan.fmax # center freq of top channel

    f1 = fch1
    f2 = f1 + fdmt_band
    #assert f2 < plan.fmax, f'F2 is past fmax f2={f2} fmax={fmax}' # not sure this is ncessary
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
    ocff = fdmt.cff(f1, fmin, fmax, fmin) + 0.5/float(plan.dmax) # Need to add a bit extra. not sure why.
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
        expectshape = (plan.nuvrest, plan.nt, plan.ndout, plan.nuvwide)
        assert blk.shape == expectshape, f'Invalid input blk shape. Was {blk.shape} expected {expectshape}'
        npix = self.plan.npix
        g = np.zeros((npix, npix), dtype=self.plan.dtype)
        fdmt_band = plan.ncin*plan.foff # Bandwidth of FDMT
        sum_max = 0
        sum_predict = 0

        for irun, fdmtrun in enumerate(plan.fdmt_plan.fdmt_runs):
            mincell = min(fdmtrun, key=lambda cell:cell.chan_start)
            minchan = mincell.chan_start
            fch1 = mincell.freqs[0]
            blkdm = int(np.round(idm*idm_cff(fch1, plan)))
            toff = int(np.round(idm*offset_cff(fch1, plan)))
            blkt = 2*t - toff
            logging.debug('Gridder idm=%s t=%s irun=%s  minchan=%s blkdm=%s toff=%s blkt=%s fch1=%s idm_cff=%s off_cff=%s', idm, t, irun, minchan, blkdm, toff,  blkt, fch1, idm*idm_cff(fch1, plan), idm*offset_cff(fch1, plan))

            for iuv, uvcell in enumerate(fdmtrun):
                upix, vpix = uvcell.uvpix
                if (blkt >= 0):
                    v1 = blk[irun, blkt+0, blkdm, iuv]
                else:
                    v1 = complex(0,0)

                if (blkt+1 >= 0):
                    v2 = blk[irun, blkt+1, blkdm, iuv]*1j
                else:
                    v2 = complex(0,0)


                # This is summing because some UV Cells get split to do the FDMT and we need to recombine them
                g[vpix, upix] += v1 + v2

                #g[npix-1-vpix, npix-1-upix] += np.conj(v1) - np.conj(v2)
                # Assume we're gridding with DC bin at [0,0] rather than [npix/2, npix/2]
                g[-vpix,-upix] += np.conj(v1) - np.conj(v2)

        logging.debug('Gridding idm=%s t=%s sum_max=%s sum_predict=%s', idm, t, sum_max, sum_predict)
                
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
        logging.info('Saving %s candidates to %s',len(self.candidates), fname)
        fmt = '%f' if len(self.candidates) == 0 else '%d %d %d %d %d %0.1f'
        np.savetxt(fname, self.candidates, header='idm t xpix ypix boxwidth sn', fmt=fmt)

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
        assert nt % 2 == 0
        outshape = (nd, nt//2, npix, npix)
        logging.info("Input shape is %s. Writing output data to %s shape is (nd,nt/2,npix,npix)=%s", blk.shape, outfname, outshape)
        fout = open(outfname, 'w')
        gout = open(outgridname, 'w')
        
        assert nt % 2 == 0, 'Nt must be divisible by 2 as were doing complex-to-real gridding'

        dosave = self.plan.values.save
        if dosave:
            img_all = np.zeros(outshape, dtype=np.complex64)
            grid_all = np.zeros(outshape, dtype=np.complex64)
        else:
            img_all = None
            grid_all = None
        
        for idm in range(nd):
            for t in range(nt//2):
                #g = gridder(d[idm, 2*t, :], d[idm, 2*t+1, :])
                g = gridder(idm, t, blk)
                #g.tofile(gout)
                img = imager(g).astype(np.complex64)
                if img_all is not None:
                    img_all[idm, t, :, :] = img
                    grid_all[idm, t, :, :] = g
                #img.tofile(fout)
                # send real part in to the boxcar & grouping
                c1 = grouper(idm, 2*t, boxcar(idm, img.real))
                c2 = grouper(idm, 2*t + 1, boxcar(idm, img.imag)) # send imaginary part into p
                rlabel = 'real idm={} t={}'.format(idm, t)
                ilabel = 'imag idm={} t={}'.format(idm, t)
                logging.debug('img.real idm=%d t=%d %s', idm, t, craco.printstats(img.real, rlabel))
                logging.debug('img.imag idm=%d t=%d %s', idm, t, craco.printstats(img.imag, ilabel))
                logging.debug('grid.real idm=%d t=%d %s', idm, t, craco.printstats(g.real, 'grid.real'))
                logging.debug('grid.imag idm=%d t=%d %s', idm, t, craco.printstats(g.imag, 'grid.imag'))
                logging.debug('idm=%s t=%d t1 candidates=%d t2 candidates=%d', idm, t, len(c1), len(c2))
                if self.values.show_image:
                    fig, ax = pylab.subplots(2,2)
                    imshow_complex(ax[0,:], img, 'image idm={} t={}'.format(idm, t))
                    imshow_complex(ax[1,:], g, 'grid idm={} t={}'.format(idm, t))
                    pylab.show()

            logging.info('Finished idm=%s', idm)
    
        grouper.to_file(candname)
        fout.close()
        gout.close()
        return img_all, grid_all
        

def savefile(fname, arr):
    logging.info('Saving file %s shape=%s dtype=%s', fname, arr.shape, arr.dtype)
    np.save(fname, arr)

class CracoPipeline(Kernel):
    def __init__(self, values):
        uvsource = uvfits.open(values.uv)
        plan = craco_plan.PipelinePlan(uvsource, values)
        super(CracoPipeline, self).__init__(uvsource, plan, values)
        self.prepare = Prepare(self.uvsource, self.plan, self.values)
        self.fdmt = MiniFdmt(self.uvsource, self.plan, self.values)
        self.image = ImagePipeline(self.uvsource, self.plan, self.values)

    def __call__(self):
        for blkt, d in enumerate(self.uvsource.time_blocks(self.plan.nt)):
            if self.plan.values.save:
                savefile('blk_{}_input.npy'.format(blkt), craco.bl2array(d))

            dprep = self.prepare(d)

            if self.plan.values.save:
                savefile('blk_{}_prepare.npy'.format(blkt), dprep)

            blk = self.fdmt(dprep)

            if self.plan.values.save:
                savefile('blk_{}_fdmt.npy'.format(blkt), blk)

            img, grid = self.image(blk)
            if self.plan.values.save:
                savefile('blk_{}_img.npy'.format(blkt), img)
                savefile('blk_{}_grid.npy'.format(blkt), grid)

        logging.info('Pipeline finished')

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
