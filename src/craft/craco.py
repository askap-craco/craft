#!/usr/bin/env python
"""
CRACO utilities

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import warnings
from astropy.coordinates import SkyCoord, Angle
try:
    from numba import njit, prange
except:
    print('Couldnt import numba. Bahh humbug')
    # define pass-through decorators
    def njit(f):
        return f

    def prange(x):
        return range(x)
          

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def triangular(n):
    '''
    Compute the Nth triangular number including diagaonal.
    i.e. number of baselines (including autos) with N antennas
    t = N*(N+1)/2
    
    :param: n - n (>= 0)
    :returns: nth triangular number

    >>> triangular(1)
    1
    
    >>> triangular(2)
    3

    >>> triangular(6)
    21
    '''
    t = n*(n+1)//2

    return t

def triangular_index(x, y, n, raster='xy'):
    '''
    Returns the index in the triangle given x and y coordinates (x>=y) and 
    the size of the array (n)

    If raster ='xy' (the default) - it assumes we raster accross x first, then down y
    if raster = 'yx', its the opposite raster order. I.e. y first, then x

    Assumes 0,0 is the top-left of the image

    >>> triangular_index(0,0,5)
    0

    >>> triangular_index(3,0,5)
    3

    >>> triangular_index(1,1,5)
    5

    >>> triangular_index(4,3,5)
    13

    >>> triangular_index(4,4,5)
    14

    >>> triangular_index(106, 71, 256)
    15726

    >>> triangular_index(0,0,5, 'yx')
    0

    >>> triangular_index(3,0,5,'yx')
    6

    >>> triangular_index(1,1,5,'yx')
    2

    >>> triangular_index(4,3,5,'yx')
    13

    >>> triangular_index(4,4,5,'yx')
    14



    '''
    assert raster in ('xy','yx'), f'Invalid raster value:{raster}'

    assert 0 <= x < n, 'Invalid values 0<= x(%d) < n(%d)'%(x, n)
    assert 0 <= y < n , 'Invalid values 0<= y(%d) < n(%d)'%(y, n)
    assert x >= y , 'x(%d) > y(%d)'% (x, y)

    if raster == 'xy':
        i = triangular(n) - triangular(n - y) + x - y
    else:
        i = triangular(x) + y
        
    
    return i

def conjidx(i, npix):
    '''
    Returns conjugate index of i
    i.e. if index is 0, returns 0, else returns npix-1
    '''
    cidx = 0 if i == 0 else npix - i
    return cidx


def conjuv(uv, npix):
    '''
    Returns conjugate uv pixel coordinate
    '''
    cuv = ((conjidx(uv[0], npix), conjidx(uv[1], npix)))
    
    return cuv
    

def make_upper(uvpix, npix):
    '''
    Make the given uvpix tuple upper hermetian.
    i.e. always returns a tuple with u >= v
    assumes (0,0) pixel is top left and raster order is rows then columns

    TODO: FIX TESTS - as they're rubbish!
    
    >>> make_upper((2,1), 256)
    (2, 1)
    
    >>> make_upper((1,1), 256)
    (1, 1)

    >>> make_upper(255-3, 255-1), 256)
    (3, 1)
    '''
    
    u, v = uvpix
    if u >= v:
        return (u, v)
    else:
        return conjuv(uvpix, npix)


def bl2ant(bl):
    '''
    Convert baseline to antena numbers according to UV fits convention
    Antenna numbers start at 1 and:

    baseline = 256*ant1 + ant2

    :see: http://parac.eu/AIPSMEM117.pdf

    :returns: (ant1, ant2) as integers

    >>> bl2ant(256*1 + 2)
    (1, 2)

    >> bl2ant(256*7 + 12)
    (7, 12)
    '''
    ibl = int(bl)
    a1 = ibl // 256
    a2 = ibl % 256

    assert a1 >= 1
    assert a2 >= 1

    return (a1, a2)

def ant2bl(ants):
    '''
    Convert antenna tuple into UV baseline index using UV fits convention
    Antenna numbers start at 1.

    :param: ants is tuple (a1, a2) where a1 and a2 are 1-based antenna numbers
    :see: bl2ant

    >>> ant2bl((7,12)) == 256*7 + 12
    True
    '''
    a1, a2 = ants
    assert a1 >= 1
    assert a2 >= 1
    ibl = a1*256 + a2
    return ibl


def bl2array(baselines):
    '''
    Converts baseline dictionary into an array sorted by baseline id
    :returns: np array of shape [nbl, nf, nt]
    If input is a masked array, returns a masked array too
    '''
    blids = sorted(baselines.keys())
    nbl = len(blids)
    indata = baselines[blids[0]]
    tfshape = indata.shape
    fullshape = [nbl]
    fullshape.extend(tfshape)

    d = np.zeros(fullshape, dtype=np.complex64)

    if np.ma.isMaskedArray(indata):
        d = np.ma.masked_array(d, mask=np.zeros(fullshape, dtype=bool), fill_value=indata.fill_value)
    
    for idx, blid in enumerate(blids):
        d[idx, ...] = baselines[blid]

    return d

def get_bl_length(baselines, blid):
    '''
    Gets the baseline distance for a given baseline
    :params: baselines (dict), blid (float)
    :returns: Distance in seconds units (float)
    '''
    dist = np.sqrt(baselines[blid]['UU']**2 + baselines[blid]['VV']**2 + baselines[blid]['WW']**2)
    return dist

def runidxs(x):
    ''' 
    Return the indexes of the start an end of a list numbers that might be equal

    '''
    istart = 0
    #for i in xrange(1, len(x)):
    for i in range(1, len(x)):
        if x[i] != x[istart]:
            yield (istart, i-1)
            istart = i
            
    yield (istart, i)

def arcsec2rad(strarcsec):
    return np.radians(float(strarcsec)/3600.)

def get_max_uv(baselines, fmax):
    '''
    Get maximum values of u,v in lambda
    '''
    ulam_max = max([abs(bldata['UU'])*fmax for bldata in list(baselines.values())])
    vlam_max = max([abs(bldata['VV'])*fmax for bldata in list(baselines.values())])

    assert ulam_max > 0 and vlam_max > 0, f'Shouldnt have 0 uv. {ulam_max} {vlam_max} {fmax} {baselines}'
    
    return (ulam_max, vlam_max)


def image_fft(g, scale='none'):
    '''
    Do the complex-to-complex imaging FFT with the correct shifts and correct inverse thingy
    If g.shape = (Npix, Npix) then we assume the center of the UV plane is at
    Npix/2, Npix/2 = DC
    Noramlised by np.prod(img.shape) - which I think is the same as the Xilinx FFT
    
    :scale: 'none' or None for raw FFT output. 'prod' for np.prod(g.shape)

    '''
    # The old version was incorrect!
    #cimg = np.fft.fftshift(np.fft.ifft2(g)).astype(np.complex64)

    if scale == 'none':
        s = 1.0
    elif scale == 'prod':
        s = np.prod(g.shape)
    
    cimg = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g)))/s
    return cimg

def printstats(d, prefix=''):
    '''
    Print and return statistics
    :prefix: prefix to stats
    :d: np.array to find stats of
    :returns: string describing statistics
    '''
    maxidx = d.argmax()
    maxpos = np.unravel_index(maxidx, d.shape)
    dmax = d.max()
    dmin = d.min()
    dmean = d.mean()
    dstd = d.std()
    dsum = d.sum()
    sn = 0 if dstd == 0 else dmax/dstd
    s = '{prefix} max/min/mean/rms/sum/S/N = {dmax:.2e}/{dmin:0.2e}/{dmean:0.2e}/{dstd:0.2e}/{dsum:0.2e}/{sn:0.1f} peak at {maxpos}'.format(**locals())
    #print s
    return s

def coord2lm(coord:SkyCoord, phase_center:SkyCoord):
    '''
    Convert given skyCoord to LM suitable for pointsource()
    :returns: lenght 2 numpy array of direction cosines (l,m)
    '''
    theta = coord.dec.rad - phase_center.dec.rad
    psi = np.cos(coord.dec.rad)*(phase_center.ra.rad - coord.ra.rad)
    lm = np.sin([psi, theta])
    return lm

def psitheta2coord(psitheta, phase_center:SkyCoord):
    '''
    Compute coordinate of source given offsets in psi/theta and a phase center as SkyCoord
    '''
    psi,theta = psitheta
    expected_dec = phase_center.dec + theta 

    # RA needs to be decremented by source cos dec
    expected_ra = np.degrees(phase_center.ra.rad - psi.rad/np.cos(expected_dec.rad))
    expected_pos = SkyCoord(expected_ra, expected_dec, unit='deg')
    return expected_pos

def pointsource(amp, lm, freqs, baseline_order, baselines, noiseamp=0):
    '''
    Returns simulted visibilities for a point source with given amplitude at given value of lm =(l, m) in radians
    offset from the phase center
        
    :amp: amplitude
    :lm: tuple of l,m as direction cosines i.e. l=sin(psi), m=sin(theta)
    :freqs: list of frequencies (Hz)
    :baseline_order: array of NBL of desired order of baselinse
    :baselines: dictionary keyed by blid containign UVWs
    :noiseamp: add nois with tgiven ampliutude
    :returns: np.array of complex dtype with shape (nbl, nchan)
    '''
    nbl = len(baselines)
    nf = len(freqs)
    assert np.all(freqs > 500e6), 'Invalid frequencies'
    
    l, m = lm
    dout = np.empty((nbl, nf), dtype=np.complex64)
    for ibl, blid in enumerate(baseline_order):       
        # baselines in seconds
        uvw_sec = np.array(baselines[blid][:3])
        
        # convert UVW coordinates to wavelegths
        u = uvw_sec[0]*freqs
        v = uvw_sec[1]*freqs
        w = uvw_sec[2]*freqs

        # TMS equation 3.7 - don't include 1/sqrt(1 - l*l - m*m) denomiator term for point sources#
        # TMS has a minus sign
        # Andre Offringa's simulator has a plus sign
        # go figure?
        # Without the minus sign also matches miriad
        vis = amp*np.exp(2j*np.pi*(u*l + v*m + w*(np.sqrt(1.0 - l*l - m*m) - 1.0)))
        if noiseamp > 0:
            vishape = vis.shape
            noise = noiseamp*(np.random.randn(*vishape) + 1j*np.random.randn(*vishape))
            vis += noise

        dout[ibl, :] = vis

    return dout


def time_blocks_with_uvws(vis, nt, flagant=[], flag_autos=True, mask=False, fetch_uvws = True):
    '''
    Generator that returns nt time blocks of the given visiblity table

    :flagant: list of antenna numbers. If an antenna number is in this list, all baselines to this antenna 
    are removed from the output
    :flag_autos: If true, autos will be removed from the output
    :mask: If true, returns np.masked_array with the masked values where the weights in the file were 0
    :fetch_uvws: If True, extracts the UVWs for each time bin and returns them as a dict. If False, returns an empty dict
    :returns: Tuple of 1)Dictionary- keyed by baseline ID, value = np.array size (nchan, nt) dtype=complex64,
                       2)Dictionary - keyed by baseline ID, value = np.array size (3, nt), dtype=float64 (empty if fetch_uvws is False)
    :see: bl2ant()
    '''

    nrows = vis.size
    inshape = vis[0].data.shape
    nchan = inshape[-3]
    npol = inshape[-2]

    shape = (nchan, npol, nt)
    uvw_shape = (3, nt)
    logging.info('returning blocks for nrows=%s rows nt=%s visshape=%s', nrows, nt, vis[0].data.shape)
    d = {}
    uvws = {}
    t = 0
    d0 = vis[0]['DATE']
    first_blid = None
    #for irow in xrange(nrows):
    for irow in range(nrows):
        row = vis[irow]
        blid = row['BASELINE']

        a1,a2 = bl2ant(blid)
        if a1 in flagant or a2 in flagant or (flag_autos and a1 == a2):
            continue

        if first_blid is None:
            first_blid = vis[0]['BASELINE']
            first_blid_again = False
        else:
            first_blid_again = blid == first_blid
            
        date_changed = row['DATE'] > d0
        
        #logging.(irow, blid, bl2ant(blid), row['DATE'], d0, t)
        if first_blid_again or date_changed: # integration finifhsed when we see first blid again. date doesnt have enough resolution
            t += 1
            tdiff = row['DATE'] - d0
            d0 = row['DATE']
            logging.debug('Time change or baseline change irow=%d, len(d)=%d t=%d d0=%s rowdate=%s tdiff=%0.2f millisec. First bl again? %s date changed? %s', irow, len(d), t, d0, row['DATE'], tdiff*86400*1e3, first_blid_again, date_changed)

            if t == nt:
                logging.debug('Returning block irow=%d, len(d)=%d t=%d d0=%s rowdate=%s tdiff=%0.2f millisec', irow, len(d), t, d0, row['DATE'], tdiff*86400*1e3)
                yield d, uvws
                d = {}
                uvws = {}
                t = 0


        if blid not in list(d.keys()):
            dvalue = np.zeros(shape, dtype=np.complex64)
            if fetch_uvws:
                uvw_value = np.zeros(uvw_shape, dtype=np.float64)
                uvws[blid] = uvw_value

            if mask:
                dvalue = np.ma.masked_array(dvalue, mask=np.zeros(shape, dtype=bool), copy=False, fill_value=0+0j)
            
            d[blid] = dvalue

        db = d[blid]

        if mask:
            # mask values if input row weight is zero
            db[..., t].data.real = row.data[...,0]
            db[..., t].data.imag = row.data[...,1]
            db[..., t].mask = row.data[...,2] == 0
        else:
            db[..., t].real = row.data[...,0]
            db[..., t].imag = row.data[...,1]

        if fetch_uvws:
            uvws[blid][..., t] = row['UU'], row['VV'], row['WW']
        
    if len(d) > 0:
        if t < nt - 1:
            warnings.warn('Final integration only contained t={} of nt={} samples len(d)={} nrows={}'.format(t, nt, len(d), nrows))
            # Returns without yielding - this is the end
        else:
            assert t == nt -1
            yield d, uvws

def _vis2nbl(vis):
    """
    functions to infer the number of the baselines based on the visibility data from uvfits.
    This don't consider any flagging etc as we don't care about that
    """
    d0 = vis[0]["DATE"]

    nbl = 0
    for visrow in vis:
        d = visrow["DATE"]
        if d > d0: return nbl
        nbl += 1


def time_block_with_uvw_range(
        vis, trange, flagant=[], flag_autos=True, mask=False, fetch_uvws=True,
    ):
    """
    function to fetch a range of vis data based on the trange

    Params
    ----------
    vis: astropy.io.fits.hdu.groups.GroupData (craco.uvfits.vis)
        visibility data loaded from the uvfits data
    trange: tuple => (int, int)
        starting sample index and ending sample index
    flagant: list
        a list of antenna to be flagged
    flag_autos: bool, True by default
        returned time block does not contain auto correlation if True
    mask: bool, default False
        mask out zero values if True 
    fetch_uvws: bool, default True
        fetch uvw values from the vis file if True

    Returns
    ----------
    d: numpy.ndarray or numpy.ma.core.MaskedArray
        visilibity data extracted from the vis 
    uvw: numpy.ndarray
        uvw values for the corresponding time range
    (_sstart, _ssend): two integers, tuple


    Notes
    ----------
    1) We don't care about any __nstart in this case, 
    the indices in the timestamp are all referring to the index in the vis data 
    2) this run as fast as the generator with next function, I will not remove this...
    """
    ### we need to figure out how many baselines are actually there in vis data...
    nbl = _vis2nbl(vis)
    tt = vis.size // nbl # timesample in total

    ### perform some initial checking on trange...
    _sstart, _send = trange

    if _sstart < 0:
        warnings.warn(f"wrong starting sample index - {_sstart}, change it to zero...")
        _sstart = 0
    if _send > tt:
        warnings.warn(f"wrong ending sample index - {_send}, change it to the final integration {tt}")
        _send = tt

    ### only consider the data after this time stamp
    nrows = vis.size
    inshape = vis[0].data.shape
    nchan = inshape[-3]
    npol = inshape[-2]
    nt = _send - _sstart + 1 # include both starting and ending

    shape = (nchan, npol, nt)
    uvw_shape = (3, nt)
    logging.info('returning blocks for nrows=%s rows nt=%s visshape=%s', nrows, nt, vis[0].data.shape)
    d = {}
    uvws = {}
    t = 0
    d0 = vis[0]['DATE']
    first_blid = None
    for irow in range(_sstart * nbl, nrows):
        row = vis[irow]
        blid = row['BASELINE']

        a1,a2 = bl2ant(blid)
        if a1 in flagant or a2 in flagant or (flag_autos and a1 == a2):
            continue

        if first_blid is None:
            first_blid = vis[0]['BASELINE']
            first_blid_again = False
        else:
            first_blid_again = blid == first_blid
            
        date_changed = row['DATE'] > d0
        
        #logging.(irow, blid, bl2ant(blid), row['DATE'], d0, t)
        if first_blid_again or date_changed: # integration finifhsed when we see first blid again. date doesnt have enough resolution
            t += 1
            tdiff = row['DATE'] - d0
            d0 = row['DATE']
            # if t is greater than nt, jump out of the loop
            if t >= nt: return d, uvws, (_sstart, _send)
            logging.debug('Time change or baseline change irow=%d, len(d)=%d t=%d d0=%s rowdate=%s tdiff=%0.2f millisec. First bl again? %s date changed? %s', irow, len(d), t, d0, row['DATE'], tdiff*86400*1e3, first_blid_again, date_changed)
        if blid not in list(d.keys()):
            dvalue = np.zeros(shape, dtype=np.complex64)
            if fetch_uvws:
                uvw_value = np.zeros(uvw_shape, dtype=np.float64)
                uvws[blid] = uvw_value

            if mask:
                dvalue = np.ma.masked_array(dvalue, mask=np.zeros(shape, dtype=bool), copy=False, fill_value=0+0j)
            
            d[blid] = dvalue

        db = d[blid]

        if mask:
            # mask values if input row weight is zero
            db[..., t].data.real = row.data[...,0]
            db[..., t].data.imag = row.data[...,1]
            db[..., t].mask = row.data[...,2] == 0
        else:
            db[..., t].real = row.data[...,0]
            db[..., t].imag = row.data[...,1]

        if fetch_uvws:
            uvws[blid][..., t] = row['UU'], row['VV'], row['WW']

    return d, uvws, (_sstart, _send)

def time_blocks(vis, nt, flagant=[], flag_autos=True, mask=False):
    d_uvw = time_blocks_with_uvws(vis, nt, flagant, flag_autos, mask, fetch_uvws = False)
    for d in d_uvw:
        yield d[0]

def grid(uvcells, Npix):
    '''
    Grid the data onto an Npix grid
    '''
    np2 = int(float(Npix)/2.)
    g = np.zeros((Npix, Npix), dtype=np.complex)

    for b in uvcells:
        upix, vpix = b.uvpix
        g[vpix, upix] += b.nchan
        g[Npix-vpix, Npix-upix] += b.nchan

    return g

def fdmt_transpose(dblk, ncu=1, nt_per_cu=2):
    '''
    Transpose the given data in to a format suitable for the Image kernel
    
    :dblk: Data block. Shape = (nuv, ndm, nt) - complex valued. nt >= ncu*nt_per_cu
    :ncu: Number of processing Compute Units
    :nt_per_cu: Number of times samples processed in parallel by a single FFT compute unit. 
    Usually this = 2 as a CU grids 2 time samples onto a complex plane to produce 2 FFT outputs (real/imag)
    :returns: data reshaped to (nuv, ndm, nt/(ncu*nt_per_cu), ncu, nt_per_cu)
    '''
    
    nuv, ndm, nt = dblk.shape
    assert ncu >= 1
    assert nt_per_cu >= 1
    assert np.iscomplexobj(dblk)
    assert nt >= ncu*nt_per_cu, 'Invalid tranpose'
    
    # Add CU axis
    # // Input order is assumed to be [DM-TIME-UV][DM-TIME-UV]
    #// Each [DM-TIME-UV] block has all DM and all UV, but only half TIME
    #// in1 is attached to the first block and in2 is attached to the second block
    #// The first half TIME has [1,2, 5,6, 9,10 ...] timestamps
    #// The second half TIME has [3,4, 7,8, 11,12 ...] timestamps
    nt_rest = nt//(ncu*nt_per_cu)

    # this is NUV, NDM, NT/(NCU*NTPERCU), NCU, NTPERCU order
    rblk = dblk.reshape(nuv, ndm, nt_rest, ncu, nt_per_cu)

    # Tranpose to (NCU, NDM, NT/(NCU*NTPERCU), NTPERCU, NUV) order
    outorder = (3, 1, 2, 4, 0)
    oblk = np.transpose(rblk, outorder)

    assert oblk.shape == (ncu, ndm, nt_rest, nt_per_cu, nuv), 'Invalid shape = {}'.format(oblk.shape)
    
    return oblk

def fdmt_transpose_inv(oblk, ncu=1, nt_per_cu=2, nuv=None, ndm=None, nt=None):
    '''
    Transpose the given data from a format suitable for the image kernel back into something sane.
    
    :oblk: Output block - otuput by fdmt_tarnspose() or a flattened array. If flattend, then ndm and nuv must be specified
    :ncu: Number of processing Compute Units
    :nt_per_cu: Number of times samples processed in parallel by a single FFT compute unit. 
    Usually this = 2 as a CU will grid 2 time samples onto a complex plane to produce 2 FFT outputs (real/imag)
    
    :returns: Data in sane ordering: (nuv, ndm, nt)
    '''
    assert ncu >= 1
    assert nt_per_cu >= 1
    assert np.iscomplexobj(oblk)

    if oblk.ndim == 1:
        nt_rest = nt // (ncu * nt_per_cu)
        oblk = oblk.reshape(ncu, ndm, nt_rest, nt_per_cu, nuv)
    else:
        (ncu_d, ndm, nt_rest, nt_per_cu_d, nuv) = oblk.shape
        # Check shape agrees with arguments
        assert ncu == ncu_d
        assert nt_per_cu_d == nt_per_cu

    nt = nt_rest * ncu * nt_per_cu
    assert nt == ncu*nt_per_cu*nt_rest, 'Invalid tranpose'

    assert oblk.shape == (ncu, ndm, nt_rest, nt_per_cu, nuv), 'Invalid shape = {}'.format(oblk.shape)

    # Reorder back to sane ordering - aiming for NUV, NDM, NT
    order = (4, 1, 2, 0, 3)
    rblk = np.transpose(oblk, order)
    assert rblk.shape == (nuv, ndm, nt_rest, ncu, nt_per_cu), 'Invalid shape={}'.format(rblk.shape)
    
    dblk = rblk.reshape(nuv, ndm, nt)
    assert dblk.shape == (nuv, ndm, nt)
    
    return dblk

def get_freqs(hdul):
    '''
    Returns a numpy array of channel frequencies in Hz from a UVFITS HDU list

    :returns: Np array length NCHAN
    '''
    
    hdr = hdul[0].header
    fch1 = hdr['CRVAL4']
    foff = hdr['CDELT4']
    ch1 = hdr['CRPIX4']
    #assert ch1 == 1.0, f'Unexpected initial frequency: {ch1}'
    assert foff > 0, 'cant handle negative frequencies anywhere athe moment foff=%f' % foff
    vis = hdul[0].data
    nchan = vis[0].data.shape[-3]

    # Need to add 1 due to stupid FITS convention. Grr.
    freqs = (np.arange(nchan, dtype=np.float) - ch1 + 1)*foff + fch1 # Hz

    return freqs

class BaselineCell(object):
    '''
    Tells you which baseline is going to end up in a uv cell

    :blid: Baslined ID - see blid2ant
    :uvpix: (2 tuple) of (u, v) integer pixels
    :chan_start, chan_end: First and last channels inclusive
    :freqs: array of frequencies
    :npix: Number of pixels
    '''
    def __init__(self, blid, uvpix, chan_start, chan_end, freqs, npix):
        self.blid = blid
        self.uvpix = uvpix
        self.chan_start = chan_start
        self.chan_end = chan_end
        assert self.chan_end >= self.chan_start
        self.freqs = freqs
        assert len(freqs) == self.chan_end - self.chan_start + 1, 'Invalid frequency mapping channels=%d-%d nfreq=%d'% (chan_start, chan_end, len(freqs))
        self.a1, self.a2 = bl2ant(blid)
        self.npix = npix

    @property
    def chan_slice(self):
        return slice(self.chan_start, self.chan_end+1)

    @property
    def fch1(self):
        return self.freqs[0]

    @property
    def uvpix_upper(self):
        '''
        Returns the uv pixel coordinates tuple guaranteed to be in the 
        upper half plane
        
        If the supplied uvpix is in the lower half, then the (u, v)
        values will be swaped, and 'is_lower()' will return True
        
        :returns: (u, v) where u >= v always.
        '''
            
        if self.is_upper:
            retuv = self.uvpix
        else:
            retuv = conjuv(self.uvpix, self.npix)

        return retuv

    @property
    def uvpix_lower(self):
        '''
        Returns uV pixel coordinates tuple guaranteed to be in the lower half plane
        '''
        return conjuv(self.uvpix_upper, self.npix)
    

    @property
    def upper_idx(self):
        '''
        Returns the triangular upper index of this guy
        '''
        u, v = self.uvpix_upper
        i =  triangular_index(u, v, self.npix)
        
        return i

    @property
    def lower_idx(self):
        '''
        Returns the triangular index of the lower coordinates
        '''
        u, v = self.uvpix_upper
        i = triangular_index(v, u, self.npix)
        
        return i

    @property
    def is_lower(self):
        '''
        Returns True if the uvpix coordinates supplied in the constructor
        where in the lower half plane. I.e. if u < v
        '''
        return not self.is_upper

    @property
    def is_upper(self):
        '''
        Returns True if the uvpix coordinates supplied in the constructor
        where in the upper half plane including diagonal. I.e. if u >= v
        '''
        u, v = self.uvpix
        return u >= v

    @property
    def nchan(self):
        return self.chan_end - self.chan_start + 1

    @property
    def chan_slice(self):
        '''
        Returnns a slice that you can use to slice out the channels of interest
        i.e. start:end+1
        '''
        return slice(self.chan_start, self.chan_end+1)

    def extract(self, baseline_block):
        cstart = self.chan_start
        cend = self.chan_end+1
        # pad with zeros
        alld = baseline_block[self.blid]
        padded_d = np.zeros_like(alld)
        padded_d[cstart:cend, :] = alld[cstart:cend, :]
        return padded_d

    def __str__(self):
        s = 'Cell blid=%s chan=%d-%d freq=%f-%f uvpix=%s upper_idx=%s uvpix_upper=%s' % (self.blid, self.chan_start, self.chan_end, self.freqs[0], self.freqs[-1], self.uvpix, self.upper_idx, self.uvpix_upper)
        return s

    __repr__ = __str__

def baseline2uv(plan, baseline_data, uv_data):
    '''
    Create UV data from baseline data given the plan
    Asumes uv_data is zeroed appropriately.

    :plan: PipelinePlan from fdmt_plan
    :baseline_data: dict keyed by baselineId, each containing (nchan, nt) nparray - see time_blocks()
    :uv_data: np array shape (plan.nuvrest, plan.nt, plan.ncin, plan.nuvwide)

    i.e. if plan is new or changed, set uv_data to zero. Otherwise you can re-use it.

    '''
    
    uv_shape = (plan.nuvrest, plan.nt, plan.ncin, plan.nuvwide)
    assert uv_data.shape == uv_shape, f'Invalid uv_data shape. Was {uv_data.shape} expected {uv_shape}'

    for irun, run in enumerate(plan.fdmt_plan.runs):
        for iuv, uv in enumerate(run.cells):
            #print(f'run {irun}/{len(plan.fdmt_plan.runs)}')
            blid = uv.blid
            cstart = uv.chan_start
            cend = uv.chan_end+1
            nchan = cend - cstart # number of channels being copied from this baseline
            run_cstart = run.chan_start # first channel in this run
            out_cstart = cstart - run_cstart # The channel in this run that the data will be copied into
            out_cend = out_cstart + nchan
            u = iuv
            urest = irun
            bldata = baseline_data[blid] # get the baseline

            #print(blid, u, urest, cstart, cend,  nchan, run_cstart, out_cstart, out_cend)
            
            assert 0 <= out_cstart <= plan.ncin
            assert out_cstart <= out_cend <= plan.ncin
            assert bldata.shape == (plan.nf, plan.nt)

            # copy all nt from the baseline at the give cstart:cend channel range
            # into out_cstart:out_cend in the uv_data
            uv_data[urest, :, out_cstart:out_cend, u] = bldata[cstart:cend, :].T

    return uv_data

@njit
def baseline2uv_numba(lut, baseline_data, uv_data):
    nrun, nuvwide, _ = lut.shape
    for irun in prange(nrun):
        for iuv in prange(nuvwide):
            blidx, cstart, cend, out_cstart, out_cend, do_conj = lut[irun, iuv, :]
            if blidx == -1:
                break

            if do_conj:
                uv_data[irun, :, out_cstart:out_cend, iuv] = np.conj(baseline_data[blidx, cstart:cend, :].T)
            else:
                uv_data[irun, :, out_cstart:out_cend, iuv] = baseline_data[blidx, cstart:cend, :].T

class FastBaseline2Uv:
    def __init__(self, plan, conjugate_lower_uvs=False):
        '''
        Numba-compiled version of baseline2uv assuming data has been smashed with bl2array  pre-compiled indexes

        :plan:PipelinePLan to operate on
        :conjugate_lower_uvs: If True, conjugate lower UV data
        '''
        self.plan = plan
        # initialise with -1 - if those values are -1 in the execution code, then we quite the loop
        self.lut = np.ones((len(plan.fdmt_plan.runs), plan.nuvwide, 6), np.int16)*-1
        blids = sorted(plan.baselines.keys())

        for irun, run in enumerate(plan.fdmt_plan.runs):
            for iuv, uv in enumerate(run.cells):
                #print(f'run {irun}/{len(plan.fdmt_plan.runs)}')
                blid = uv.blid
                cstart = uv.chan_start
                cend = uv.chan_end+1
                nchan = cend - cstart # number of channels being copied from this baseline
                run_cstart = run.chan_start # first channel in this run
                out_cstart = cstart - run_cstart # The channel in this run that the data will be copied into
                out_cend = out_cstart + nchan
                u = iuv
                urest = irun
                blidx = blids.index(blid) 
            
                assert 0 <= out_cstart <= plan.ncin
                assert out_cstart <= out_cend <= plan.ncin
                #bldata = plan.baselines[blid]
                #assert bldata.shape == (plan.nf, plan.nt)
                #uv_data[urest, :, out_cstart:out_cend, u] = bldata[cstart:cend, :].T
                do_conj = int(conjugate_lower_uvs and uv.is_lower)
                self.lut[irun, iuv, :] = [blidx, cstart, cend, out_cstart, out_cend, do_conj]

    def __call__(self, baseline_data, uv_data):
        '''
        Convert baselines to UV data

        :baseline_data: baseline data sorted into an array: see bl2array. Shape=(nbl, nc, nt)
        :uv_data: output uv data shape  (nurest, nt, ncin, nuvwide)
        '''

        assert uv_data.shape == self.plan.uv_shape, f'Invalid uv_data shape. Was {uv_data.shape} expected {self.uv_shape}'
        assert baseline_data.shape == self.plan.baseline_shape, f'Invalid baseline_data shape. Was {baseline_data.shape} expected {self.plan.baseline_shape}'

        baseline2uv_numba(self.lut, baseline_data, uv_data)


    
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
