#!/usr/bin/env python
"""
Cleans multibeam filterbanks using eigenflagging/subtraction.

See Kocz 2012 and Wax and Kailath 1985

Copyright (C) CSIRO 2018
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from . import sigproc
import scipy.stats
import warnings
from . import dada

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

last_mean = None
last_std = None

def flagplot(ax, df, thresh, values, title):
    ax.plot(df[values.plot_beam, 0, :])
    if not np.isinf(thresh):
        ax.axhline(thresh, c='r')
        ax.axhline(-thresh, c='r')

    ax.set_title(title)

def flag(d, values):
    '''
    Flag data similarly to fredda
    '''
    global last_mean
    global last_std
    nbeam, nint, nchan = d.shape

    dmean = d.mean(axis=1)[:, np.newaxis, :]
    dstd = d.std(axis=1)[:, np.newaxis, :]

    if last_mean is None:
        last_mean = dmean
                 
    if last_std is None:
        last_std = dstd
        
    meanoff = dmean/last_mean - 1
    stdoff = dstd/last_std - 1

    mflags = abs(meanoff) > values.mean_thresh
    sflags = abs(stdoff) > values.std_thresh

    dkurt = scipy.stats.kurtosis(d, axis=1)[:, np.newaxis, :]
    kflags = abs(dkurt) > values.kurt_thresh

    #gtest fabs(mean*mean/variance/nsamps_per_int - rescale_dtype(1));
    gtest = abs(dmean*dmean/dstd/dstd/float(values.nsamps_per_int) - 1)
    gflags = gtest > values.gtest_thresh

    allflags = ~(mflags | sflags | kflags | gflags)

    dout = d * allflags
    #dout[allflags] = 0
    last_mean = dmean
    last_std = dstd

    if values.show:
        fig, ax = pylab.subplots(2,3)
        ax = ax.flatten()

        flagplot(ax[0], meanoff, values.mean_thresh, values, 'meanoff')
        flagplot(ax[1], stdoff, values.std_thresh, values, 'stdoff')
        flagplot(ax[2], dkurt, values.kurt_thresh, values, 'kurt')
        flagplot(ax[3], gtest, values.gtest_thresh, values, 'gtest')
        ax[4].imshow(dout[values.plot_beam, :, :].T, aspect='auto')


    return dout

def eigenflag(d, values):
    assert values.in_rescale, 'You must rescale the input before you eigenflag, otherwise youll get very rubbish results'
    nbeam, nint, nchan = d.shape
    dout = np.empty_like(d)
    vic = np.empty((nchan, nbeam))
    uic = np.empty((nchan, nbeam, nbeam))
    k = values.nsub
    for ic in range(nchan):
        r = d[:, :, ic]
        thecov = np.cov(r)
        #dcovar[ic, :, :] = thecov
        try:
            u, v, vh = np.linalg.svd(thecov)
            vic[ic, :] = v
            uic[ic, :, :] = u
            for t in range(nint):
                rfi = np.dot(u.T, r[:, t])
                #rfiout[:, t, ic] = rfi[:]
                rfi[:k] = 0
                dout[:, t, ic] = np.dot(rfi, u.T)
        except ValueError: # Happens when svd fails
            warnings.warn('SVD Failed chan={}'.format(ic))
            dout[:,:,ic] = d[:,:,ic]


    if values.show:
        fig, ax = pylab.subplots(2,2)
        ax = ax.flatten()
        ax[0].plot(vic)
        ax[0].set_ylabel('Singular value')
        ax[0].set_xlabel('Channel number')
        ax[1].imshow(uic[:, :, 0].T, aspect='auto')
        ax[1].set_ylabel('Beam')
        ax[1].set_xlabel('Channel number')
        ax[1].set_title('Principal eigenvector')
        ax[2].plot(d[:, :, values.plot_channel].T)
        ax[3].plot(dout[:,:,values.plot_channel].T)

        pylab.show()


    return dout

def subtract_dm0(d, values):
    dm0 = d.mean(axis=2)
    dout = d - dm0[:,:,np.newaxis]

    return dout

def eigenflag_dm0(d, values):
    nbeam, nint, nchan = d.shape
    dout = np.empty_like(d)
    dm0 = d.mean(axis=2)
    thecov = np.cov(dm0)
    u, v, vh = np.linalg.svd(thecov)
    k = values.nsub
    opt_dm0 = np.empty((nbeam, nint))

    for t in range(nint):
        r = dm0[:, t]
        rfi = np.dot(u.T, r)
        rfi[k:] = 0
        rfixed = np.dot(rfi, u.T)
        opt_dm0[:, t] = rfixed

    dout = d - opt_dm0[:,:,np.newaxis]
        
    if values.show:
        # find most correlated beams - peak off diagonal element
        maxidx = np.argmax(abs(np.triu(thecov, k=1)))
        maxbeams = np.unravel_index(maxidx, thecov.shape)
        fig, ax = pylab.subplots(2,5)
        fig.suptitle('DM0 subtraction largest beam correlation {}'.format(maxbeams))
        #ax[0].imshow(dm0.T, aspect='auto')
        #ax[1].imshow(opt_dm0.T, aspect='auto')
        #ax[0,0].plot(dm0.T)
        #ax[0,1].plot(opt_dm0.T)
        ax[0,0].plot(dm0.T[:, maxbeams])
        ax[0,0].set_title('DM0 before cleaning')

        newdm0 = dout.mean(axis=2)
        ax[0,1].plot(newdm0.T[:, maxbeams])
        print('Have I subtracted everything', dm0.T[:, maxbeams].max(), newdm0.T[:, maxbeams].max())
        ax[0,1].set_title('DM0 after cleaning')

        ax[0,2].plot(v)
        ax[0,3].imshow(thecov)
        ax[1,0].imshow(d[values.plot_beam,:,:].T, aspect='auto')
        ax[1,1].imshow(dout[values.plot_beam,:,:].T, aspect='auto')
        ax[1,2].imshow(u)
        ax[1,3].plot(u[:,:3])# first eigenvector

        #ax[0,4].scatter(dm0[maxbeams[0],:], dm0[maxbeams[1],:])

        ax[0,4].plot(opt_dm0.T[:, maxbeams])
        ax[0,4].set_title('The DM0 RFI')

        ax[1,4].plot(dm0.std(axis=1), 'o')
        ax[1,4].plot(opt_dm0.std(axis=1), 'o')
        pylab.show()

    return dout

def rescale(d, values):
    nbeam, nint, nchan = d.shape
    m = d.mean(axis=1)[:, np.newaxis, :]
    s = d.std(axis=1)[:, np.newaxis, :]
    dout = np.zeros_like(d)
    dout = (d - m)/s
    mask =~np.isfinite(dout)
    dout[mask] = 0
    
    return dout

class SigprocReader(object):
    def __init__(self, files, values):
        self.insfs = [sigproc.SigprocFile(f) for f in values.files]
        self.nsamp_total = min([s.nsamples for s in self.insfs])
        self.nsamp = values.nsamp
        self.nblocks = self.nsamp // values.nsamp
        startmjd = np.array([s.tstart for s in self.insfs])
        assert np.all(startmjd == startmjd[0])
        self.nbeams = len(self.insfs)
        self.values = values

    def header(self, beamid):
        return self.insfs[beamid].header

    def outfile(self, beamid):
        return self.insfs[beamid].filename

    def read(self, blkid):
        start = blkid*self.nsamp
        end = start  + self.nsamp
        d = np.array([s[start:end] for s in self.insfs])
        return d



class DadaReader(object):
    def __init__(self, fin, values):
        self.infile = dada.DadaFile(fin)
        self.nblocks = self.infile.nblocks
        self.nbeams = int(self.infile.hdr.get_value('NBEAM'))
        self.values = values

    def header(self, beamid):
        ''' Returns a dictionary suitable for use as filterbank header'''
        h = self.infile.hdr
        hdr = {}
        bra = list(map(float, h.get_value('BEAM_RA').split(',')))[beamid]
        bdec = list(map(float, h.get_value('BEAM_DEC').split(',')))[beamid]
        hdr['foff'] = float(h.get_value('BW'))
        hdr['fch1'] = float(h.get_value('FREQ'))
        hdr['tstart'] = float(h.get_value('MJD_START'))
        hdr['tsamp'] = float(h.get_value('TSAMP'))
        hdr['nifs'] = 1
        hdr['nchans'] = int(h.get_value('NCHAN'))
        hdr['nbits'] = self.infile.dtype.itemsize*8
        hdr['src_raj'] = bra
        hdr['src_dej'] = bdec

        return hdr

    def outfile(self, beamid):
        return 'beam.{:02d}.fil'.format(beamid)

    def read(self, blkid):
        d = self.infile.get_block(blkid)
        
        # discard first dimension
        # input dimension is [1,nbeam,nchan,nint]
        # output dimesion should be [nbeam, nint, nchan]
        d = d[0, :, :, :].transpose(0, 2, 1)
        
        return d

        
def open_files(values):
    r = None
    if len(values.files) > 1 and all([s.endswith('.fil') for s in values.files]):
        r = SigprocReader(values.files, values)
    elif len(values.files) == 1 and values.files[0].endswith('.dada'):
        r = DadaReader(values.files[0], values)
    else:
        raise ValueError('Unknown file formats: {}'.format(values.files))

    return r

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-n','--nsamp', type=int, help='Number of samples per block to calculate covariance over and subtract', default=256)
    parser.add_argument('-k', '--nsub', type=int, default=1, help='Fixed Number of eigenvalues to subtract (todo: Dynamic)')
    parser.add_argument('-D','--outdir', required=True, help='Output directory for cleaned filterbanks')
    parser.add_argument('-f', '--flag', action='store_true', help='Do zero flagging')
    parser.add_argument('-e','--eigenflag', action='store_true', help='Apply eigenflagging')
    parser.add_argument('-s','--subtract-dm0', action='store_true', help='Subtract dm0')
    parser.add_argument('-d','--eigenflag-dm0', action='store_true', help='Eigenflag on dm0')
    parser.add_argument('-r','--in-rescale', action='store_true', help='Rescale input')
    parser.add_argument('--start-samp', type=int, help='Sample to start on', default=0)
    parser.add_argument('-N', '--nblock', type=int, help='Number of blocks to process')
    parser.add_argument('-M', '--mean-thresh', type=float, help='Relative mean threshold', default=np.inf)
    parser.add_argument('-T', '--std-thresh', type=float, help='Relative std threshold', default=np.inf)
    parser.add_argument('-K', '--kurt-thresh', type=float, help='Kurtosis threshold', default=np.inf)
    parser.add_argument('-G', '--gtest-thresh', type=float, help='Gtest threshold', default=np.inf)
    parser.add_argument('-I','--nsamps-per-int', type=int, help='Samples per integration for gtest threshold', default=1024)
    
    parser.add_argument('--show', action='store_true', help='Show plots')
    parser.add_argument('-b','--plot-beam', help='Beam index to plot', type=int, default=0)
    parser.add_argument('-c','--plot-channel', help='Channel index to plot', type=int, default=0)
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    reader = open_files(values)
    nbeams = reader.nbeams
    outfile_names = [os.path.join(values.outdir, reader.outfile(b)) for b  in range(nbeams)]
    for f in outfile_names:
        d = os.path.dirname(f)
        if not os.path.isdir(d):
            os.makedirs(d)

    outsfs = [sigproc.SigprocFile(f, mode='w', header=reader.header(b)) for b,f in enumerate(outfile_names)]
    nsamp = values.nsamp
    startblock = values.start_samp // nsamp
    nblock = values.nblock

    if values.nblock is None:
        endblock = reader.nblock
    else:
        endblock = startblock + values.nblock
    
    k = values.nsub

    print(startblock, endblock)

    for b in range(startblock, endblock):
        d = reader.read(b)
        dtype = d.dtype

        d = d.astype(float)

        if values.flag:
            d = flag(d, values)

        if dtype == np.uint8:
            d -= 128
            d /= 18

        if values.in_rescale:
            d = rescale(d, values)

        if values.eigenflag:
            d = eigenflag(d, values)
 
        if values.eigenflag_dm0:
            d = eigenflag_dm0(d, values)

        if values.subtract_dm0:
            d = subtract_dm0(d, values)

        if dtype == np.uint8:
            d *= 18
            d += 128

        for iout, sout in enumerate(outsfs):
            dout = d[iout, :, :].astype(dtype).flatten()
            if isinstance(dout, np.ma.MaskedArray):
                dout = dout.data
            
            dout.tofile(sout.fin)


    for s in outsfs:
        s.fin.close()

if __name__ == '__main__':
    _main()
