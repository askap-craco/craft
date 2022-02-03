#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys
import numpy as np
import pylab
import pickle as pickle


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Dumps interesting infromation about a craft_fil.py pickle dump')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('files', help='Files to process', nargs='+')
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

        
    for fin in values.files:
        #dump_str(fin, values)
        plot(fin, values)

def plot(fin_name, values):
    fin = open(fin_name, 'rb')
    hdr = pickle.load(fin)
    for k in sorted(hdr.keys()):
        print(k, hdr[k])

    ints = []
    xs= []
    ys = []
    fAs = []
    fBs = []
    bAs = []
    bBs = []

    nignore = 10

    while True:
        try:
            (i, x, y, fA, fB, bA, bB) = pickle.load(fin)
            fA = np.array(fA)
            print(x.shape)

            if np.all(fA==0):
                print('Zero packet i', i)
                continue

            if nignore > 0:
                nignore -= 1
                continue


            ints.append(i)
            xs.append(x)
            ys.append(y)
            fAs.append(fA)
            fBs.append(fB)
            bAs.append(bA)
            bBs.append(bB)
        except EOFError as UnpicklingError:
            break


    ints, xs, ys, fAs, fBs, bAs, bBs = list(map(np.array, (ints, xs, ys, fAs, fBs, bAs, bBs)))

    print('All ints', set(ints.flatten()))

    #Arrays are in time, beam, freq order

    print(xs.shape, ys.shape, ints.flatten().shape)
#    fig, (beamax, intax) =  pylab.subplots(2,1)
#    beamax.plot(beams)
#    beamax.set_ylabel('Beam')
#    intax.plot(ints)
#    intax.set_ylabel('Int')
#    intax.set_xlabel('packet number')

    pylab.figure()
    pylab.plot(ys[:,0,:], 'x') # Frequencies each beam are always identical
    pylab.title('Frequencies')
    pylab.xlabel('packet number')

    pylab.figure()
    pylab.plot(10*np.log10(xs[:,0,:]))
    pylab.title('Powers (dB)')
    pylab.xlabel('Packet number')

    fig, (fa_ax, fb_ax, fd_ax) = pylab.subplots(3,1)
    fa_ax.plot(fAs[:,0,:])
    fa_ax.set_ylabel('Start frame number')
    fb_ax.plot(fBs[:,0,:])
    fb_ax.set_ylabel('Stop frane number')
    fd_ax.plot(fBs[:,0,:] - fAs[:,0,:])
    fd_ax.set_ylabel('Frame number difference')

    # beam 0 int 0 
    #ntime, nbeams, nfreqs = x.shape
    # why is nints 3?
    pylab.figure()
    pylab.plot((fAs[1:,0,:] - fAs[0:-1,0,:])/int(hdr['CRAFT_INT_TIME']))
    pylab.title('Frane number increment - beam0 All freqs')
    pylab.show()


def dump_str(fin_name, values):
    fin = open(fin_name, 'rb')
    hdr = pickle.load(fin)
    for k in sorted(hdr.keys()):
        print(k, hdr[k])

        
    while True:
        try:
            (beam, i, x, y, fA, fB, bA, bB) = pickle.load(fin)
            print('b',beam, 'i',i, 'len',len(x))
            print('x', x)
            print('fA', fA)
            print('fB', fB)
            print('bA', bA)
            print('bB', bB)
        except EOFError as UnpicklingError:
            break

    

def bat2mjd(bat):
    utcdt = askap.time.bat2utcDt(bat)
    mjd = askap.time.utcDt2mjd(utcdt)

if __name__ == '__main__':
    _main()
