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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"
fname_re = re.compile('.*ak(.*).*.(\d\d).fil*')

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

    for d in values.files:
        plot(d, values)

def plot(d, values):
    
    all_files = []
    
    all_data = []
    for dirpath, dirnames, filenames in os.walk(d):
        for f in filenames:
            if f.endswith('.cand'):
                fname = os.path.join(dirpath, f)
                if os.path.getsize(fname) != 0:
                    all_files.append(fname)
                    cand = np.loadtxt(fname)
                    if cand.ndim == 1:
                        cand.shape = (1, len(cand))

                    all_data.extend(cand)

    all_data = np.array(all_data)
    cand = all_data

    if len(cand) == 0:
        print d, 'is empty'
        return 


    fig, axes = pylab.subplots(2,2)
    fax = axes.flatten()

    print cand.shape

    snr = cand[:, 0]
    dm = cand[:, 5]
    t = cand[: ,2]
    width = cand[:, 3]
    mask = (snr > 10.) 

    snr = cand[mask, 0]
    dm = cand[mask, 5]
    t = cand[mask, 2]
    width = cand[mask, 3]

    fig.text(0.5, 0.98, 'Heimall {}'.format(d), ha='center', va='top')

    fax[0].scatter(t, dm)
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

    fout = 'heim_{}.png'.format(d).replace(r'/','')
    print 'Saving', fout
    fig.savefig(fout)
    

    if values.show:
        pylab.show()


        
    pylab.close()
            

                                 


        
    

if __name__ == '__main__':
    _main()
