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

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

# Indeces of the center, min and max columns for time and dm
snidx =0
t0 = 1
t1 = 7
t2 = 8
d0 = 4
d1 = 9
d2 = 10
count = 11

header = 'S/N, sampno, secs from file start, boxcar, idt, dm, beamno, sampno_start, sampno_end, idt_start, idt_end, ncands'
intf = '%d'
floatf = '%0.3f'
formats = (floatf, intf, floatf, intf, intf, floatf, intf, intf, intf, intf, intf, intf)

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-t', '--tdist', type=int, help='Sample distance', default=32)
    parser.add_argument('-d','--ddist', type=int, help='Idt distance', default=20)
    parser.add_argument('-p', '--plot', action='store_true', help='Show plots', default=False)
    
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)



    for fin in values.files:
        d = np.loadtxt(fin)
        if len(d) == 0:
            logging.info('%s is empty', fin)
            continue
        assert d.shape[1] == 7
        # make new array. Concatenate time twice and dm twice at the end, along with a counter
        hstack = (d, d[:, t0:t0+1], d[:, t0:t0+1], d[:, d0:d0+1], d[:, d0:d0+1], np.ones((d.shape[0], 1)))
        dnew = np.hstack(hstack)
        dout = fof(dnew, values)
        np.savetxt(fin+'.fof', dout, fmt=formats, header=header)
        logging.info('%s reduced from %d to %d candidates', fin, dnew.shape[0], dout.shape[0])


def errplot(alldarr, **kwargs):
    alldarr = np.array(alldarr)
    xerr = np.array([-alldarr[:, t1] + alldarr[:, t0], alldarr[:, t2] - alldarr[:, t0]])
    yerr = np.array([-alldarr[:, d1] + alldarr[:, d0] , alldarr[:, d2] - alldarr[:, d0]])
    return pylab.errorbar(alldarr[:, t0], alldarr[:, d0], xerr=xerr, yerr=yerr, ls='none', **kwargs)


def find_friends_mask(myd, d, values):
    mysamp1, mysamp2, myidt1, myidt2 = myd[t1], myd[t2], myd[d1], myd[d2]
    mysamp = myd[t0]
    myidt = myd[d0]
    assert mysamp1 <= mysamp <= mysamp2
    assert myidt1 <= myidt <= myidt2
    tmask1 = mysamp1 - values.tdist <= d[:, t2]
    tmask2 = mysamp2 + values.tdist >= d[:, t1]

    dmask1 = myidt1 - values.ddist <= d[:, d2]
    dmask2 = myidt2 + values.ddist >= d[:, d1]
    
    tmask = tmask1 & tmask2
    dmask = dmask1 & dmask2
    fullmask= tmask & dmask
    return fullmask

def fof(d, values):
    alld = []
    last_numd = d.shape[0]
    old_d = d.copy()
    iteration = 1
    while True:
        bestd, drest = fof_iter(d, values)
        alld.append(bestd)

        d = drest
        if drest.shape[0] == 0:
            if values.plot:
                pylab.figure()
                errplot(alld, color='r', marker='o')
                errplot(old_d, marker='x')
                pylab.title('Iteration %d' % iteration)
                pylab.savefig('Iteration%d.png'%iteration)

            iteration += 1

            d = np.array(alld)
            if d.shape[0] == last_numd:
                break

            last_numd = d.shape[0]
            old_d = d.copy()
            alld = []

    pylab.show()

    return d
        


def fof_iter(d, values):
    myd = d[0, :]
    fullmask = find_friends_mask(myd, d, values)
    related = d[fullmask, :]

    assert related.shape[0] != 0

    if related.shape[0] == 1: # we're the only thing nearby
        best_d = myd
    else: # update points we ate
        # find the highest sn
        best_idx = np.argmax(related[:, snidx])
        best_d = related[best_idx, :]
        
        best_d[t1] = np.min(related[:, t1])
        best_d[t2] = np.max(related[:, t2])
        
        best_d[d1] = np.min(related[:, d1])
        best_d[d2] = np.max(related[:, d2])
        best_d[count] += related.shape[0]

        assert best_d[t1] <= best_d[t0] <= best_d[t2] 
        assert best_d[d1] <= best_d[d0] <= best_d[d2]



    return best_d, d[~fullmask, :]


    
    

if __name__ == '__main__':
    _main()