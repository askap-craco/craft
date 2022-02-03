#!/usr/bin/env python
"""
Sets up craft and downlaods channel mapping

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys
import askap.craft.beamformer 
from askap.craft.beamformer import CraftBeamformer, FilterTypes, ZoomModes
from askap.craft.crafthdr import DadaHeader
from askap.craft.freqconfig import FreqConfig
from askap.craft.cmdline import strrange
import pylab
import numpy as np

def _main():
    from argparse import ArgumentParser
    from askap.craft.cmdline import strrange
    parser = ArgumentParser(description='Sets up craft and downloads channel mapping')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='file', help='Header file to plot')
    parser.add_argument('-c','--cards', help='List of cards')
    parser.add_argument('--show', action='store_true', help='Show plot')

    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    h = DadaHeader.fromfile(values.file)
    if values.cards:
        cards = strrange(values.cards)
    else:
        cards = None

    fconfig = FreqConfig.load_from_dada_header(h, cards)
    nbf, nchans = fconfig.freqs.shape

    print(fconfig)

    print('NUM_BEAMFORMERS {}'.format(nbf))

    assert cards is None or nbf == len(cards)

    print('FREQ {}'.format(fconfig.freq))
    print('BW {}'.format(fconfig.bw))
    print('NCHAN {}'.format(fconfig.nchan))

    fig, axes = pylab.subplots(2,1)
    ax1, ax2 = axes.flatten()


    for i in range(nbf):
        freqs = fconfig.freqs[i, :]
        chanmap = fconfig.chanmaps[i, :]
        freqmap = fconfig.freqmaps[i, :]

        if cards is None:
            cardno = i+1
        else:
            cardno = cards[i]

        print('BEAMFORMER{}_CARDNO {}'.format(i, cardno))
        print('BEAMFORMER{}_CHANMAP {}'.format(i,  ','.join(map(str, chanmap))))
        print('BEAMFORMER{}_FREQMAP {}'.format(i,  ','.join(map(str, freqmap))))
        print('BEAMFORMER{}_FREQS {}'.format(i,  ','.join(map(str, freqs))))

        if values.show:
            ax1.plot(freqs, 'o', label='ibf {}'.format(i))
            ax2.plot(freqs,  np.ones(len(freqs))*i, 'x', label='Card %d' % (i+1))



    if values.show:
        pylab.legend()
        pylab.xlabel('Channel number')
        pylab.ylabel('Frequency (MHz)')
        pylab.savefig(values.file+'.png')
        pylab.show()
            

if __name__ == '__main__':
    _main()
