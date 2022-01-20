#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
from astropy.io import fits
from .uvfits2fil import get_antname
import shutil
import pylab

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


class UvFitsReader(object):
    def __init__(self, fin):
        hdu = fits.open(fin, mode='update')
        g = hdu['PRIMARY']
        hdu.info()
        hdr = g.header
        self.fcent = hdr['CRVAL4'] / 1e6 # MHz
        self.foff = hdr['CDELT4'] / 1e6 # MHz
        self.ref_channel = hdr['CRPIX4'] # pixels
        self.g = g
        self.nchan = self.g.data[0]['DATA'].shape[3]
        self.freqs = (np.arange(self.nchan) - self.ref_channel + 1)*self.foff + self.fcent # not sure if this is out by half a channel
        self.fch1 = self.freqs[0]
        self.hdu = hdu

    def mjds(self):
        jds = self.g.data['DATE'].astype(np.float64)
        if '_DATE' in self.g.columns.names:
            dayfrac = self.g.data['_DATE'].astype(np.float64)
            fulljd = jds + dayfrac
        else:
            fulljd = jds
        
        mjd = fulljd - 2400000.5
        return mjd
        
    def __len__(self):
        return len(self.g.data)

    def __getitem__(self, i):
        row = self.g.data[i]
        bl = int(row['BASELINE'])
        ia1 = bl % 256 - 1
        ia2 = bl // 256 - 1
        a1 = get_antname(self.hdu, ia1)
        a2 = get_antname(self.hdu, ia2)
        data = row['DATA']
        nchan = data.shape[3]
        spec = np.zeros(nchan, np.complex64)
        spec.real = data[0,0,0,:,0,0]
        spec.imag = data[0,0,0,:,0,1]
        weights = data[0,0,0,:,0,2]
        jd = row['DATE']
        inttime = row['INTTIM']
        dayfrac = row['_DATE']
        fulljd = float(jd) + float(dayfrac)
        mjd = fulljd - 2400000.5

        return (a1, a2, mjd, spec, weights)

    def set_weights(self, i, weights):
        self.g.data[i]['DATA'][0,0,0,:,0,2] = weights.astype(np.float32)

    def flush(self):
        self.hdu.flush()

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='Show plot')
    parser.add_argument('-o','--outfile', help='Output file', required=True)
    parser.add_argument('-c','--candfile', help='Candidate file')
    parser.add_argument('-w','--width', help='Gaussian width represented by boxcar=0 (milliseconds)', type=float, default=1.0)
    parser.add_argument('--update', action='store_true', help='Multply FRB weight to existing weights', default=False)
    parser.add_argument('--weight-data', action='store_true', help='Weight the data')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    candidates = np.loadtxt(values.candfile)
    if candidates.ndim == 1:
        candidates.shape = (1, -1)

    print(candidates.shape)

    if candidates.shape[0] > 1:
        raise ValueError('Can\'t handle more than one candidate')
    
    cand = candidates[0, :]
    fin = values.files[0]
    fout = values.outfile
    if not os.path.exists(fout):
        logging.info('Copying to output %s', fout)
        shutil.copy(fin, fout)
        
    reader = UvFitsReader(fout)
    freqs = reader.freqs/1e3 # GHz
    print('FREQS max/min/mean,0,-1', freqs.max(), freqs.min(), freqs.mean(), freqs[0], freqs[-1])
    cdm = cand[5]
    cmjd = cand[7]
    cboxcar = cand[3]
    dm_delay_ms_max = 4.15*cdm*(freqs.max()**-2 - 0)

    cand_mjds = cmjd + 4.15*cdm*(freqs**-2 - 0)/1e3/86400. # assme ref freq infinity
    cand_mjdstart = cand_mjds.min()
    cand_mjdend = cand_mjds.max()
    mjds = reader.mjds()
    width_ms = (cboxcar+1.0) *values.width/2.0
    width_days = width_ms/1e3/86400.
    
    logging.info('CANDIDATE MJDstart %s %s %s %s %s', cmjd, cand_mjdstart, cand_mjdend, dm_delay_ms_max, freqs.max())
    logging.info('Candidates is from {}-{} seconds into file'.format((cand_mjdstart - mjds[0])*86400, (cand_mjdend - mjds[0])*86400.))

    assert mjds[0] == mjds.min()
    assert mjds[-1] == mjds.max()
    
    logging.info('UVstats mjdstart {} {}'.format(mjds.min(), mjds.max()))
    logging.info('FREQS [0]={} freqs[-1]={}'.format(freqs[0], freqs[-1]))
    row_range = np.where((mjds >= cand_mjdstart - width_days*3) & (mjds <= cand_mjdend + width_days*3))[0]


    row_start = row_range.min()
    row_end = row_range.max()
    row_range = row_end - row_start

    logging.info('ROWs %s %s %s', row_start, row_end, (row_end - row_start))
    
    total_update = 0
    pulse_mjds = mjds[row_start:row_end]
    

    weights = np.zeros((len(pulse_mjds), len(freqs)))
    for irow, mjd in enumerate(pulse_mjds):
        dm_delay_ms = 4.15*cdm*(freqs**-2 - 0) # assume reference frequency of 0
        tdiff_ms = (mjd - cmjd)*86400.*1e3
        toff = (tdiff_ms - dm_delay_ms)
        #print irow, mjd, tdiff_ms, dm_delay_ms.min(), dm_delay_ms.max(), toff.min(), toff.max(), width_ms

        frb_weights = np.exp(-toff**2/(width_ms**2))

        #idxs = abs(toff) < width_ms
        #frb_weights = 1
        frb_weights[frb_weights<1e-3] = 0
        weights[irow, :] = frb_weights


    if values.weight_data:
        reader.g.data['DATA'][row_start:row_end,0,0,0,:,0,0] *= weights
        reader.g.data['DATA'][row_start:row_end,0,0,0,:,0,1] *= weights
        reader.g.data['DATA'][row_start:row_end,0,0,0,:,0,2] = 1
        reader.g.data['DATA'][0:row_start, 0,0,0,:,0,:] =0
        reader.g.data['DATA'][row_end:, 0,0,0,:,0,:] = 0

    else:
        reader.g.data['DATA'][0:row_start, 0,0,0,:,0,2] = -1
        reader.g.data['DATA'][row_end:, 0,0,0,:,0,2] = -1
        reader.g.data['DATA'][row_start:row_end,0,0,0,:,0,2] = weights
        

    reader.flush()

    if values.show:
        pylab.imshow(weights.T, aspect='auto')
        pylab.show()


def apply_weight2():
        
    for irow in range(row_start, row_end):
        a1, a2, mjd, spec, weights = reader[irow]
        if irow % 1000 == 0:
            print('Row {} of {} = {:0.1f}% Updated {}'.format(irow, len(reader), float(irow - row_start)/float(row_range)*100., total_update))

        tdiff_ms = (cmjd - mjd)*86400.*1e3
        frb_weights = np.exp(-(tdiff_ms - dm_delay_ms)**2/(2*c**2))
        frb_weights[frb_weights < 1e-3] = 0
        nupdate = sum(frb_weights != 0)
        total_update += nupdate

        if values.update:
            new_weights = frb_weights*weights
        else:
            new_weights = frb_weights
        
        if nupdate > 0:
            reader.set_weights(irow, new_weights)
            logging.debug('Updated {}-{} {} new weights. Tdiff_ms {} dm_delay_ms'.format(a1,a2, nupdate, tdiff_ms, dm_delay_ms))


if __name__ == '__main__':
    _main()
