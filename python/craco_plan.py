#!/usr/bin/env python
"""
Plans a CRACO scan

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import craco
import uvfits
import warnings

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def get_uvcells(baselines, uvcell, freqs, Npix):
    uvcells = []

    ucell, vcell = uvcell
    for blid, bldata in baselines.iteritems():
        #UU, VV WW are in seconds
        ulam = bldata['UU'] * freqs
        vlam = bldata['VV'] * freqs
        pix_offset = float(Npix)/2.0
        
        upix = np.round(ulam/ucell + pix_offset).astype(int)
        vpix = np.round(vlam/vcell + pix_offset).astype(int)
        if np.any((upix < 0) | (upix >= Npix) | (vpix < 0) | (vpix >= Npix)):
            warnings.warn('Pixel coordinates out of range')

        #if values.show:
        #    pylab.plot(ulam/1e3, vlam/1e3)

        uvpos = list(zip(upix, vpix))
        for istart, iend in craco.runidxs(uvpos):
            uvpix = uvpos[istart]
            assert uvpos[istart] == uvpos[iend]
            b = craco.BaselineCell(blid, uvpix, istart, iend, freqs[istart:iend+1], Npix)
            uvcells.append(b)

    uvcells = sorted(uvcells, key=lambda b:b.upper_idx)
    return uvcells

def calc_overlap(blcell, minchan, ncin):
    assert blcell.chan_start >= minchan, 'invalid minimum chan {} {}'.format(blcell.chan_start, minchan)
    maxchan = minchan + ncin
    overlap = maxchan - blcell.chan_start
    return overlap

def split_cells(cells, minchan, ncin):
    '''
    Given a list of UV cells and a channel range given by minchan and ncin, produces 2 lists
    The first list contains Cells that fint entirely within the chanel range. The second list
    contains cells outside the channel range

    :cells: List of UV cells
    :minchan: first channel in channel range
    :ncin: number of channels in range
    :returns: Tuple of lists
    '''
    full_cells = []
    leftover_cells = []
    endchan = minchan + ncin - 1 # final channel, inclusive

    for c in cells:
        assert c.chan_start >= minchan
        if c.chan_start >= minchan and c.chan_end <= endchan:
            full_cells.append(c)
        else:
            # split the cell in two parts
            npix = c.npix
            logging.debug('c %d-%d freq=%f-%f endchan=%d', c.chan_start, c.chan_end, c.freqs[0], c.freqs[-1], endchan)
            endidx = endchan - c.chan_start
            assert 0 <= endidx < ncin
            c1 = craco.BaselineCell(c.blid, c.uvpix, c.chan_start, endchan, c.freqs[:endidx+1], npix)
            c2 = craco.BaselineCell(c.blid, c.uvpix, endchan+1, c.chan_end, c.freqs[endidx+1:], npix)
            full_cells.append(c1)
            leftover_cells.append(c2)

    return full_cells, leftover_cells
    

class FdmtPlan(object):
    def __init__(self, uvcells, nuvwide=8, ncin=32):
        uvcells_remaining = set(uvcells)# copy array
        fdmt_runs = []
        run_chan_starts = []
        run_fch1 = []
        while len(uvcells_remaining) > 0:
            logging.debug('Got %d/%d uvcells remaining', len(uvcells_remaining), len(uvcells))
            minchan = min(uvcells_remaining, key=lambda uv:(uv.chan_start, uv.blid)).chan_start
            possible_cells = filter(lambda uv:calc_overlap(uv, minchan, ncin) > 0, uvcells_remaining)
            best_cells = sorted(possible_cells, key=lambda uv:calc_overlap(uv, minchan, ncin), reverse=True)
            used_cells = best_cells[0:min(nuvwide, len(best_cells))]
            full_cells, leftover_cells = split_cells(used_cells, minchan, ncin)
            mincell = min(full_cells, key=lambda cell:cell.chan_start)
            
            run_chan_starts.append(mincell.chan_start)
            run_fch1.append(mincell.fch1)
            fdmt_runs.append(full_cells)

            total_overlap = sum([calc_overlap(uv, minchan, ncin) for uv in full_cells])
            logging.debug('minchan=%d npossible=%d used=%d full=%d leftover=%d total_overlap=%d', minchan, len(possible_cells), len(used_cells), len(full_cells), len(leftover_cells), total_overlap)
            uvcells_remaining -= set(used_cells)
            uvcells_remaining.update(leftover_cells)
            
        nruns = len(fdmt_runs)
        nuvtotal = nruns*nuvwide
        logging.info('FDMT plan has ntotal=%d of %d runs with packing efficiency %f. requires efficiency of > %f', nuvtotal, nruns, float(len(uvcells))/float(nuvtotal), float(nuvtotal)/8192.0)
        self.fdmt_runs = fdmt_runs
        self.run_chan_starts = run_chan_starts
        self.run_fch1 = run_fch1

        self.nruns = nruns
        self.nuvtotal = nuvtotal
        self.total_nuvcells = sum([len(p) for p in fdmt_runs])


class PipelinePlan(object):
    def __init__(self, f, values):
        self.values = values
        logging.info('making Plan values=%s', values)

        umax, vmax = f.get_max_uv()
        lres, mres = 1./umax, 1./vmax
        baselines = f.baselines
        nbl = len(baselines)
        freqs = f.channel_frequencies

        # Cant handle inverted bands - this is assumed all over the code. It's boring
        assert freqs.min() == freqs[0]
        assert freqs.max() == freqs[-1]
        Npix = values.npix

        if values.cell is not None:
            lcell, mcell = map(craco.arcsec2rad, values.cell.split(','))
            los, mos = lres/lcell, mres/mcell
        else:
            los, mos = map(float, values.os.split(','))
            lcell = lres/los
            mcell = mres/mos
            
        lfov = lcell*Npix
        mfov = mcell*Npix
        ucell, vcell = 1./lfov, 1./mfov
        fmax = freqs.max()
        foff = freqs[1] - freqs[0]
        lambdamin = 3e8/fmax
        umax_km = umax*lambdamin/1e3
        vmax_km = vmax*lambdamin/1e3
        
        logging.info('Nbl=%d Fch1=%f foff=%f nchan=%d lambdamin=%f uvmax=%s max baseline=%s resolution=%sarcsec uvcell=%s arcsec uvcell= %s lambda FoV=%s deg oversampled=%s',
                 nbl, freqs[0], foff, len(freqs), lambdamin, (umax, vmax), (umax_km, vmax_km), np.degrees([lres, mres])*3600, np.degrees([lcell, mcell])*3600., (ucell, vcell), np.degrees([lfov, mfov]), (los, mos))

        uvcells = get_uvcells(baselines, (ucell, vcell), freqs, Npix)
        logging.info('Got Ncells=%d uvcells', len(uvcells))
        d = np.array([(v.a1, v.a2, v.uvpix[0], v.uvpix[1], v.chan_start, v.chan_end) for v in uvcells], dtype=np.uint32)
        np.savetxt(values.uv+'.uvgrid.txt', d, fmt='%d',  header='ant1, ant2, u(pix), v(pix), chan1, chan2')
        self.fdmt_plan = FdmtPlan(uvcells, values.nuvwide, values.ncin)
        self.nd = values.ndm
        self.nt = values.nt
        self.freqs = freqs
        self.npix = Npix
        self.nbox = values.nbox
        self.boxcar_weight = values.boxcar_weight
        self.nuvwide = values.nuvwide
        self.nuvmax = values.nuvmax
        assert self.nuvmax % self.nuvwide == 0
        self.nuvrest = self.nuvmax / self.nuvwide
        self.ncin = values.ncin
        self.ndout = values.ndout
        self.foff = foff
        self.dtype = np.complex64 # working data type
        self.threshold = values.threshold
        self.nbl = nbl
        self.fdmt_scale = self.values.fdmt_scale
        self.fft_scale = self.values.fft_scale
        assert self.threshold >= 0, 'Invalid threshold'
        if (self.fdmt_plan.nuvtotal >= values.nuvmax):
            raise ValueError("Too many UVCELLS")

    @property
    def nf(self):
        '''Returns number of frequency channels'''
        return len(self.freqs)


    @property
    def fmin(self):
        '''
        Returns maximum frequency
        '''
        return self.freqs[0]

    @property
    def fmax(self):
        '''
        Returns minimum frequency
        '''
        return self.freqs[-1]
        

def add_arguments(parser):
    '''
    Add planning arguments
    '''
    parser.add_argument('--uv', help='Load antenna UVW coordinates from this UV file')
    parser.add_argument('--npix', help='Number of pixels in image', type=int, default=256)
    parser.add_argument('--os', help='Number of pixels per beam', default='2.1,2.1')
    parser.add_argument('--cell', help='Image cell size (arcsec). Overrides --os')
    parser.add_argument('--nt', help='Number of times per block', type=int, default=16)
    parser.add_argument('--ndm', help='Number of DM trials', type=int, default=16)
    parser.add_argument('--nbox', help='Number of boxcar trials', type=int, default=4)
    parser.add_argument('--boxcar-weight', help='Boxcar weighting type', choices=('sum','avg','sqrt'), default='sum')
    parser.add_argument('--nuvwide', help='Number of UV processed in parallel', type=int, default=8)
    parser.add_argument('--nuvmax', help='Maximum number of UV allowed.', type=int, default=8192)
    parser.add_argument('--ncin', help='Numer of channels for sub fdmt', type=int, default=32)
    parser.add_argument('--ndout', help='Number of DM for sub fdmt', type=int, default=32)
    parser.add_argument('--threshold', type=float, help='Threshold for candidate grouper', default=10)
    parser.add_argument('--fdmt_scale', type=float, help='Scale FDMT output by this amount', default=1.0)
    parser.add_argument('--fft_scale', type=float, help='Scale FFT output by this amount. If both scales are 1, the output equals the value of frb_amp for crauvfrbsim.py', default=10.0)
    parser.add_argument('--show-image', action='store_true', help='Show image plots')
    parser.add_argument('--show-fdmt', action='store_true', help='Show FDMT plots')
    parser.add_argument('--save', action='store_true',  help='Save data as .npy for input, FDMT and image pipeline')


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Plans a CRACO scan', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='Show plots')
    add_arguments(parser)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Loading UV coordinates from file %s ', values.uv)
    f = uvfits.open(values.uv)
    plan = PipelinePlan(f, values)

    if values.show:
        f.plot_baselines()
        pylab.figure()
        nchans = [len(bc.freqs) for bc in uvcells]
        pylab.hist(nchans, bins=50)
        pylab.xlabel('Number of channels per cell')
        pylab.ylabel('Number of cells')

        pylab.show()




        

if __name__ == '__main__':
    _main()
