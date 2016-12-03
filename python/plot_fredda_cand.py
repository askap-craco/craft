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
import plotutil
import fnmatch
import plot_allbeams
import glob

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def find_files(rootd, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(rootd):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('-f','--fname', help='Candidate filename')
    parser.add_argument('--detail', action='store_true', help='Plot detail')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False, show=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for fin in values.files:
        pylab.close()
        if os.path.isdir(fin):
            fig,axes, ncand = plot_dir(fin, values)
        else:
            fig, ax = plotutil.subplots(1,1)
            ax = pylab.gca()
            v = plot_file(fin, values, ax, title=fin)
            ncand = len(v)

        fout = fin
        if fout.endswith('/'):
            fout = fout[0:-1]
            
        fout = '%s_%s.png' % (fout, values.fname)

        if ncand > 0:
            print 'Saving', fout
            fig.savefig(fout, dpi=200)
        else:
            print 'Not saving ', fin, 'is empty'
        
        if values.show:
            pylab.show()

def plot_dir(din, values):
    candfiles = find_files(din, values.fname)
    print 'len candfiles', len(candfiles)
    
    if len(candfiles) == 0:
        fig = pylab.figure()
        fig.text(0.5, 0.5, '%s has no candidate files'%din, transform=pylab.gca().transAxes)
        return fig, pylab.gca(), []

    ncols = 4
    nrows = len(candfiles)/ncols
    if nrows*ncols < len(candfiles):
        nrows += 1

    assert nrows * ncols >= len(candfiles)
    print 'nxy', nrows, ncols
    fig, axes = plotutil.subplots(nrows, ncols, sharex=True, sharey=True)
    fig.set_size_inches([8,6])
    fig.set_dpi(300)
    axes = axes.flatten()

    ncands = 0
    for ic, f in enumerate(candfiles):
        ax=axes[ic]
        antname = f.replace(din, '').replace('C000','').replace(values.fname,'').replace('/','')
        subtitle = antname
        v = plot_file(f, values, ax, labels=False, subtitle=subtitle)
        ncands += len(v)

    fig.suptitle(din)
    fig.text(0.5, 0.05, 'Time(s)', ha='center',va='bottom')
    fig.text(0.05, 0.5, 'Delta_t', rotation=90, ha='center', va='top')

    return fig, axes, ncands


def plot_file(fin, values, ax, title=None, labels=True, subtitle=None):
    vin = np.loadtxt(fin)

    if subtitle:
        ax.text(0.05, 0.95, subtitle, ha='left', va='top', transform=ax.transAxes)

    if len(vin) == 0:
        print fin, 'is empty'
        if subtitle:
            s = 'Empty'
        else:
            s = '%s empty' % fin.replace('/', '/\n')
        ax.text(0.5, 0.5, s, ha='center', va='center', fontsize=8, transform=ax.transAxes)
        return vin

    if len(vin.shape) == 1:
        vin.shape = (1, len(vin))

    print fin, vin.shape
    ncols = vin.shape[1]
    sn = vin[:, 0]
    sampno = vin[:, 1]
    time = vin[:, 2]
    boxcar = vin[:, 3]
    dm = vin[:, 4]
    if ncols > 5:
        idm = vin[:, 5]
        beamno = vin[:, 6]
    else:
        beamno = vin[:, 5]

    ubeams = set(beamno)
    for b in sorted(ubeams):
        bmask = beamno == b
        ax.plot(time[bmask], dm[bmask]+1, marker='x', ls='None', label='Beam %d' %b, picker=3)
               
    plotutil.addpick()
    
    if labels:
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('DM (delta_t)')

    if title:
       ax.set_title(fin)

    if values.detail:
        plot_details(fin, vin, values, ax, title, labels, subtitle)

    return vin

def plot_details(fin, vin, values, ax, title, labels, subtitle):
    fdir = os.path.basename(fin)

    ncand = vin.shape[0]
    ncols = vin.shape[1]
    for c in xrange(ncand):
        cand = vin[c, :]
        plot_single_beam(fin, values, cand)
        plot_all_beams(fin, values, cand)

def plot_beams(fin, values, cand, pattern, nxy, postfix):
    sn, sampno, time, boxcar, idm, dm, beamno = cand
    
    fdir = os.path.dirname(fin)
    bpath = os.path.join(fdir, pattern)
    beamfile = glob.glob(bpath)
    tstart = int(sampno - idm - idm/2.)
    tstart = max(0, tstart)
    nsamps = int(2*idm)
    assert len(beamfile) > 0, 'Couldnt find beamfiles in path {}'.format(bpath)
    print 'got ', len(beamfile), 'beams in path', bpath, 'tstart', tstart, 'nsamps', nsamps
    p = plot_allbeams.Plotter(beamfile, nxy, tstart, nsamps)
    prefix = os.path.join(fdir, 's{:d}_b{:d}_idm{:d}_{}'.format(int(sampno), int(beamno), int(idm), postfix))
    p.draw()
    p.saveall(prefix)
    if values.show:
        pylab.show()

    p.closeall()


def plot_single_beam(fin, values, cand):
    beamno = cand[-1]
    beampat =  r'*.{:02d}.fil'.format(int(beamno))
    plot_beams(fin, values, cand, beampat, (1,1), 'singlebeam')

def plot_all_beams(fin, values, cand):
    plot_beams(fin, values, cand, '*000000.*.fil', (6,6), 'allbeams')


if __name__ == '__main__':
    _main()
