#!/usr/bin/env python
"""
Makes dedispersed visibilities suitable for the CRACO imaging pipelien that does gridding, FFT and boxcar.


Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from astropy.io import fits
from scipy import constants
import fdmt
import warnings
from craco import *

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class BaselineCell(object):
    def __init__(self, blid, uvpix, chan_start, chan_end, freqs):
        self.blid = blid
        self.uvpix = uvpix
        self.chan_start = chan_start
        self.chan_end = chan_end
        self.freqs = freqs
        self.a1, self.a2 = bl2ant(blid)

    @property
    def nchan(self):
        return self.chan_end - self.chan_start + 1

    def extract(self, baseline_block):
        cstart = self.chan_start
        cend = self.chan_end+1
        # pad with zeros
        alld = baseline_block[self.blid]
        padded_d = np.zeros_like(alld)
        padded_d[cstart:cend, :] = alld[cstart:cend, :]
        return padded_d


def fdmt_baselines(hdul, baselines, uvcells, values):
    hdr = hdul[0].header
    fch1 = hdr['CRVAL4']
    foff = hdr['CDELT4']
    ch1 = hdr['CRPIX4']
    assert ch1 == 1.0, 'Unexpected initial frequency'
    vis = hdul[0].data
    nchan = vis[0].data.shape[-3]
    freqs = np.arange(nchan)*foff + fch1 # Hz
    nrows = vis.size
    thefdmt = fdmt.Fdmt(fch1*1e6, foff*1e6, nchan, values.ndm, values.nt)
    nbl = len(baselines)
    baselines = sorted(baselines.keys())
    nuv = len(uvcells)
    for blkt, d in enumerate(time_blocks(vis, values.nt)):
        dblk = np.zeros((nuv, values.ndm, values.nt), dtype=np.complex64)
        logging.info('UV data output shape is %s blkt=%s d.shape=%s', dblk.shape, blkt, len(d))
        # FDMT everything
        for iuv, uvd in enumerate(uvcells):
            cell_data = uvd.extract(d)
            logging.debug('blkt %s iuv %s cell_data shape=%s', blkt, iuv, cell_data.shape)
            # Truncating times for the moment, as we don't keep history
            dblk[iuv, :, :] = thefdmt(cell_data)[:, :values.nt]

            if values.show:
                fig, ax = pylab.subplots(1,2)
                ax[0].imshow(abs(dblk[iuv, :, :]))
                ax[1].imshow(np.angle(dblk[iuv, :, :]))
                pylab.show()

        if values.outfile:
            outfile = values.outfile
        else:
            outfile = values.files
            
        fname = '{}.ndm{:d}_nt{:d}.b{:d}.uvdata.{}'.format(outfile, values.ndm, values.nt, blkt, values.format)
        # Transpose to [NDm, NT, NUV] order
        dblk = np.transpose(dblk, (1, 2, 0))
        logging.info('Writing shape (NDM, NT, NUV)=%s to %s in format=%s num nonozero=%d', dblk.shape, fname, values.format, np.sum(dblk != 0))
        if values.format == 'npy':
            np.save(fname, dblk)
        else:
            dblk.tofile(fname)
        

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Makes dedispersed visibilities', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--npix', help='Number of pixels in image', type=int, default=256)
    parser.add_argument('--cell', help='Image cell size (arcsec)', default='10,10')
    parser.add_argument('--nt', help='Number of times per block', type=int, default=256)
    parser.add_argument('--ndm', help='Number of DM trials', type=int, default=16)
    parser.add_argument('--outfile', help='Output filename base. Defualts to input filename')
    parser.add_argument('--format', help='Output format', choices=('npy','raw'), default='raw')
    parser.add_argument('--no-fft-shift', help='Dont shift fft', action='store_false', dest='fft_shift', default=True)
    parser.add_argument('-s','--show', help='Show plots', action='store_true')
                        
    parser.add_argument(dest='files', nargs='?')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Running %s with arguments %s', sys.argv[0], values)
    hdul = fits.open(values.files)
    vis = hdul[0].data

    # get data for first integration to work out UVW coordinatse
    d0 = vis[0]['DATE']
    baselines = {}
    for i in range(vis.size):
        row = vis[i]
        baselines[row['BASELINE']] = row
        if row['DATE'] != d0:
            break

    hdr = hdul[0].header
    assert hdr['CTYPE4'] == 'FREQ', 'Cant find freq axis'
    fch1 = hdr['CRVAL4']
    foff = hdr['CDELT4']
    ch1 = hdr['CRPIX4']
    assert ch1 == 1.0, 'Unexpected initial frequency'
    nchan = vis[0].data.shape[-3]
    freqs = np.arange(nchan)*foff + fch1 # Hz
    lcell, mcell = map(arcsec2rad, values.cell.split(','))
    Npix = values.npix
    lfov = lcell*Npix
    mfov = mcell*Npix
    ucell, vcell = 1./lfov, 1./mfov
    logging.info('Fch1=%f foff=%f uvcell=%s', fch1, foff, (ucell, vcell))
    uvcells = []

    for blid, bldata in baselines.iteritems():
        #UU, VV WW are in seconds
        ulam = bldata['UU'] * freqs
        vlam = bldata['VV'] * freqs
        pix_offset = 0
        if values.fftshift:
            pix_offset = Npix/2
        
        upix = np.round(ulam/ucell + pix_offset).astype(int)
        vpix = np.round(vlam/vcell + pix_offset).astype(int)
        if np.any((upix < 0) | (upix >= Npix) | (vpix < 0) | (vpix >= Npix)):
            warnings.warn('Pixel coordinates out of range')

        if values.show:
            pylab.plot(ulam/1e3, vlam/1e3)

        uvpos = list(zip(upix, vpix))
        for istart, iend in runidxs(uvpos):
            uvpix = uvpos[istart]
            assert uvpos[istart] == uvpos[iend]
            b = BaselineCell(blid, uvpix, istart, iend, freqs[istart:iend+1])
            uvcells.append(b)

    uvcells = sorted(uvcells, key=lambda b:b.uvpix)
    d = np.array([(f.a1, f.a2, f.uvpix[0], f.uvpix[1], f.chan_start, f.chan_end) for f in uvcells], dtype=np.uint32)
    
    np.savetxt(values.files+'.uvgrid.txt', d, fmt='%d',  header='ant1, ant2, u(pix), v(pix), chan1, chan2')
    fdmt_baselines(hdul, baselines, uvcells, values)

    if values.show:
        pylab.xlabel('U(klambda)')
        pylab.ylabel('V(klambda)')
        
        fix, ax = pylab.subplots(1,2)
        g = grid(uvcells, values.npix)
        ax[0].imshow(abs(g), aspect='auto', origin='lower')
        ax[1].imshow(image_fft(g).real, aspect='auto', origin='lower')
        pylab.show()
    
        
    
    
        

if __name__ == '__main__':
    _main()
