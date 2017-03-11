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

markers = mpl.markers.MarkerStyle().filled_markers

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('-f','--fname', help='Candidate filename')
    parser.add_argument('-o','--outfile', help='Output file name')
    parser.add_argument('-w','--max-boxcar', help='max width to plot', default=32, type=int)
    parser.add_argument('-d','--min-dm', help='minimum dm to plot', default=0, type=float)
    parser.add_argument('-b','--min-sn', help='minimum S/N to plot', default=0, type=float)

    parser.add_argument('--detail', action='store_true', help='Plot detail')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False, show=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    cmap = plt.get_cmap('rainbow')
    maxboxcar = 32.
    scalarmap = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(1., maxboxcar))
    scalarmap._A = [] # URGH http://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots


    for fin in values.files:
        pylab.close()
        if os.path.isdir(fin):
            fig, axes, ncand = plot_dir(fin, values, scmap=cmap)
        else:
            fig, ax = plotutil.subplots(1,1)
            ax = pylab.gca()
            try:
                v = plot_file(fin, values, ax, title=fin,scmap=cmap)
                ncand = len(v)
            except:
                logging.exception('Could not plot file')
                ncand = 0
            
        # add colorbar
        fig.subplots_adjust(right=0.87)
        cbar_ax = fig.add_axes([0.90, 0.1, 0.02, 0.8])
        cbar = fig.colorbar(scalarmap, cax=cbar_ax)
        cbar.set_label('Width (samples)')

        fout = fin
        if fout.endswith('/'):
            fout = fout[0:-1]
        
        if values.outfile is None:
            fout = '%s_%s.png' % (fout, values.fname)
        else:
            fout = values.outfile

        if ncand > 0:
            print 'Saving', fout
            fig.savefig(fout, dpi=200)
        else:
            print 'Not saving ', fin, 'is empty'
        
        if values.show:
            pylab.show()

def plot_dir(din, values, scmap=None):
    candfiles = find_files(din, values.fname)
    print 'len candfiles', len(candfiles)
    
    if len(candfiles) == 0:
        fig = pylab.figure()
        fig.text(0.5, 0.5, '%s has no candidate files'%din, transform=pylab.gca().transAxes)
        return fig, pylab.gca(), 0

    ncols = 4
    nrows = len(candfiles)/ncols
    if nrows*ncols < len(candfiles):
        nrows += 1

    assert nrows * ncols >= len(candfiles)
    print 'nxy', nrows, ncols
    fig, axes = plotutil.subplots(nrows, ncols, sharex=True, sharey=True)
    fig.set_size_inches([8,6])
    #fig.set_dpi(300)
    axes = axes.flatten()

    ncands = 0
    for ic, f in enumerate(candfiles):
        ax=axes[ic]
        antname = f.replace(din, '').replace('C000','').replace(values.fname,'').replace('/','')
        subtitle = antname
        v = plot_file(f, values, ax, labels=False, subtitle=subtitle, scmap=scmap)
        ncands += len(v)

    fig.suptitle(din)
    fig.text(0.5, 0.05, 'Time(s)', ha='center',va='bottom')
    fig.text(0.05, 0.5, 'DM (pc/cm3)', rotation=90, ha='center', va='top')

    return fig, axes, ncands


def plot_file(fin, values, ax, title=None, labels=True, subtitle=None, scmap=None):
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
    boxcar = vin[:, 3]
    dm = vin[:, 5]
    sn = vin[:, 0]
    mask = (boxcar < values.max_boxcar) & (dm >= values.min_dm) & (sn >= values.min_sn)
    badvin = vin[~mask, :]
    # plot failed detections in grey
    vin = vin[mask, :]
    ncols = vin.shape[1]
    sn = vin[:, 0]
    sampno = vin[:, 1]
    time = vin[:, 2]
    boxcar = vin[:, 3]
    idm = vin[:, 4]
    dm = vin[:, 5]
    beamno = vin[:, 6]
    ubeams = set(beamno)

    ax.scatter(badvin[:, 2], 1+badvin[:, 5], marker='.', alpha=0.4, edgecolors='face', c='0.5')

    for ib, b in enumerate(sorted(ubeams)):
        bmask = beamno == b
        ax.scatter(time[bmask], 1+dm[bmask], s=sn[bmask]**2, marker=markers[ib % len(markers)], c=boxcar[bmask], cmap=scmap, label='Beam %d' %b, picker=3, edgecolors='face', alpha=0.4, norm=mpl.colors.Normalize(1., 32.))
        ax.set_yscale('log')
        ax.set_ylim(1, 5000)
               
    plotutil.addpick()
    
    if labels:
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('1+DM (pc/cm3)')

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
