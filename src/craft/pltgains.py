#!/usr/bin/env python
"""
Plots gains

Copyright (C) CSIRO 2018
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import datetime
import matplotlib.dates as mdates


__author__ = "Keith Bannister <keith.bannister@csiro.au>"
firstline = ' 0 06:00:01      -3.558      9.413      0.000      5.007    -13.164    -11.763'

col1_len = len(' 0')
col2_len = len('06:00:01    ')

def readmirgain(f):
    with open(f, 'rU') as fin:
        curr_data = []
        col1_value = None
        col2_value = None
        for iline, line in enumerate(fin):
            if line.startswith('#'):
                continue
            
            col1 = line[0:col1_len].strip()
            col2 = line[col1_len:col2_len].strip()
            line_data = line[col2_len:].split()
            if col1 != '' and col2 != '': # start of a new row
                if len(curr_data) > 0: # yield old row
                    yield (col1_value, col2_value, curr_data)
                    
                col1_value = col1
                col2_value = col2
                curr_data = []

            curr_data.extend(line_data)

    if len(curr_data) > 0:
        yield (col1_value, col2_value, curr_data)
                

def load_gains(f):
    real = list(readmirgain(f))
    imag = list(readmirgain(f.replace('.real','.imag')))
    assert len(real) == len(imag)
    nrows = len(real)
    reald = []
    imagd = []
    for row in range(nrows):
        assert len(real[row][2]) == len(imag[row][2])
        reald.append(real[row][2])
        imagd.append(imag[row][2])

    g = np.array(reald).astype(float) + 1j*np.array(imagd).astype(float)
    return g


def load_file(f, values):
    g = load_gains(f)
    #ax[0].plot(abs(g), 'x')
    nsamp, nant = g.shape
    timestr,card,beam = f.split('_')[0:3]
    time = datetime.datetime.strptime(timestr, '%Y%m%d%H%M%S')

    return (time, card, beam, g)

def load_all(files, values):
    data = {}
    for f in files:
        time, card, beam, g = load_file(f, values)
        if (card, beam) not in list(data.keys()):
            data[(card,beam)] = []
            
        data[(card, beam)].append((time, g[0, :]))

    return data

    
def plot(data, values):
    lines = []
    labels = []
    size = list(map(int, values.nxy.split(',')))
    fig, ax = pylab.subplots(size[0],size[1], sharex=True, sharey=True)
    ax = ax.flatten()

    plt.locator_params(axis='x', nbins=3)
    for card, beam in list(data.keys()):
        alld = data[(card, beam)]
        all_times = np.array([d[0] for d in alld])
        all_g = np.array([d[1] for d in alld])

        nsamp, nant = all_g.shape
        for iant in range(nant):
            g = all_g[:, iant]
            gdelta = g / g[0]
            l, = ax[iant].plot(all_times, np.angle(gdelta, deg=True), marker='o', ls='-')
            ax[iant].set_title('Ant %d'%(iant+1))
            ax[iant].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            ax[iant].xaxis.set_major_locator(plt.MaxNLocator(3))
        
        lines.append(l)
        labels.append('{} {}'.format(card, beam))


    print(lines, labels)
    fig.legend(lines, labels)
    fig.text(0.01, 0.5, 'Phase change (deg)', rotation=90, va='center', ha='left')
    fig.text(0.5, 0.01, 'Time (HH:MM)', va='bottom', ha='center')


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('--nxy', help='Number of subplots x,y', default='5,6')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    data = load_all(values.files, values)
    plot(data, values)

    pylab.show()
    

if __name__ == '__main__':
    _main()
