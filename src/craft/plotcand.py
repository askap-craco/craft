#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import re
from . import craftobs
import glob
from . import heimdall


__author__ = "Keith Bannister <keith.bannister@csiro.au>"
fname_re = re.compile('.*ak(.*).*.(\d\d).fil*')

class CandidatePlotter(object):
    def __init__(self, scans):
        self._scans = scans
        self._make_plots()
        self._update(len(self._scans) - 1)

    def _clear_axes(self):
        #self._sn_ax.clf()
        #self._dm_ax.clf()
        #self._cand_ax.clf()
        #self._figmain.clf()
        # self._figbeams.clf()
        pass


    def _make_plots(self):
        self._figmain = plt.figure()
        g = (2,3)
        self._sn_ax  = plt.subplot2grid(g, (0, 0))
        self._dm_ax  = plt.subplot2grid(g, (0, 1))
        self._width_ax = plt.subplot2grid(g, (0, 2))
        self._cand_ax = plt.subplot2grid(g, (1, 0), colspan=3)
        self._figbeams, self._axbeams = pylab.subplots(6,6, sharex=True, sharey=True)

        
    def _update(self, hdrnum):

        scan = self._scans[hdrnum]
        self._scan = scan
        logging.debug("Updating hdrnum %d hdr %s", hdrnum, scan)
        cands = scan.load_candidates()
        self._cands = cands

        beams = sorted(cands.keys())
        self._clear_axes()
        for beam in beams:
            c = cands[beam]
            self._plot_cands(c)

        self._plot_single_cand(0, 0)

    def _plot_single_cand(self, ibeam, icand):
        c = self._cands[ibeam][icand, :]
        tstart = c[2]
        ntimes = 128
        beams = self._scan.load_beams(tstart, ntimes)
        for iax, ax in enumerate(self._axbeams.flat):
            lbeam = (beams[:, iax, :]).T
            ax.imshow(lbeam, aspect='auto')
        

    def _plot_cands(self, cand):
        snr = cand[:, 0]
        dm = cand[:, 5]
        t = cand[: ,2]
        width = cand[:, 3]
        mask = (snr > 10.) 
        
        snr = cand[mask, 0]
        dm = cand[mask, 5]
        t = cand[mask, 2]
        width = cand[mask, 3]
        print('plotting', len(snr), 'candidates from', cand.shape)
        if len(snr) == 0:
            return


        #fig.text(0.5, 0.98, 'Heimall {}'.format(d), ha='center', va='top')
        fax = (self._cand_ax, self._sn_ax, self._dm_ax, self._width_ax)
        fax[0].scatter(t, dm+1, marker='o')
        fax[0].grid(True)

        #fax[0].get_yaxis().set_scale('log')

        fax[0].set_xlabel('Time')
        fax[0].set_ylabel('DM')
        fax[1].hist(snr, histtype='step')
        fax[1].set_ylabel('num')
        fax[1].set_xlabel('S/N')
        fax[2].hist(dm, histtype='step')
        fax[2].set_ylabel('Count')
        fax[2].set_xlabel('DM')
        fax[3].hist(width, histtype='step')
        fax[3].set_ylabel('Count')
        fax[3].set_xlabel('width')

        

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('--show', action='store_true')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    scans = craftobs.find_scans(values.files[0])

    c = CandidatePlotter(scans)
    
    pylab.show()
    

if __name__ == '__main__':
    _main()
