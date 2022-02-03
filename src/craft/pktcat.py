#!/usr/bin/env python
"""
Convvert packet dump data from a packtdump file to a filterbank file. See craft_pktdump.py

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys

import askap.time
import numpy as np
import time
import pylab
import pickle as pickle
import datetime
import askap.craft.pktdump as pktdump
from askap.craft.sigproc import SigprocFile

        
def bat2mjd(bat):
    utcdt = askap.time.bat2utcDt(bat)
    mjd = askap.time.utcDt2mjd(utcdt)
    return mjd

def get_channel(hdr, address):
    fpga = int(np.log2(hdr['srcDevice']))
    card = int(address[0].split('.')[-1]) -1
    assert fpga >= 0 and fpga <= 5, 'Invalid FPGA %s' % fpga
    assert card >=0 and card <= 8, 'Invalid card %s' % card

    channel = card*6 + fpga
    return channel

def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='converst packet dump data to filterbank file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='file')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    fin = values.file
    f = pktdump.pktopen(fin)
    main_hdr = f.hdr
    print(main_hdr)
    intcount = main_hdr['CRAFT_INT_CYCLES']
    intime = main_hdr['CRAFT_INT_TIME']
    batint = intime*27./30.
    
    ntimes = intcount*4
    first_bat = None
    nrx = 0
    currbat = 0
    MAX_INT_COUNT = 8
    frame1 = None
    bat1 = None
    #ignore_cards = ['10.2.13.1']
    ignore_cards = []
    print_times = False
    last_times = {}

    print('CARDS', main_hdr['CARDS'])

    cards = list(map(int, main_hdr['CARDS'].split(',')))
    ncards = len(cards)
    nchans_per_fpga = 1 # change to 8 when we  get non-truncated  packets - CRAFT-6
    nfpgas_per_card = 6
    nchan = ncards*nchans_per_fpga*nfpgas_per_card
    npol = 1
    dtypes = {32: np.float32, 8: np.uint8}
    dtype = dtypes[main_hdr['nbits']]
    nslots = 8
    currslot = 0
    vdata = np.zeros((nslots, nchan, npol), dtype=dtype)
    vints = np.zeros(nslots, dtype=np.int)
    
    for hdr, address, rawdata, pktdate in f:
        if address[0] in ignore_cards:
            print('Skipping bad card ' + address[0])
            continue

        if hdr['packetType'] != 129:
            continue

        if hdr['fragmentID'] > 0:
            #print 'Skipping fragment', hdr, len(rawdata)
            continue

        times = np.fromstring(rawdata[0:8*2*MAX_INT_COUNT], dtype=np.uint64)
        data  = np.fromstring(rawdata[ntimes*8:], dtype=dtype)
        d = times
        nts = ntimes/2
        start_frames = d[0:nts:2]
        stop_frames = d[1:nts:2]
        start_bats = d[nts:2*nts:2]
        stop_bats = d[nts+1:2*nts:2]
        pktbat = start_bats[0]
        pktframe = start_frames[0]
        hdrbat = hdr['bat']
        dev = hdr['srcDevice']
        seqno = hdr['packetSequence']
        if print_times:
            print('times', times)
            print('start bats', start_bats)
            print('stop bat', stop_bats)
            print('bat diffs', stop_bats - start_bats)
            print('start frames', start_frames)
            print('stop frames', stop_frames)
            print('frame diffs', stop_frames - start_frames)

        chan = get_channel(hdr, address)
        #print start_frames - stop_frames, type(start_frames), start_frames.dtype, start_frames, (start_frames - stop_frames), (start_frames + stop_frames), start_frames[0] + stop_frames[0]
        
        fdiff = stop_frames - start_frames 
        bdiff = stop_bats - start_bats 

        if pktbat < first_bat:
            print('Skipping', first_bat, pktbat, first_bat - pktbat, hdr, len(rawdata), d)
            continue

        if frame1 is None:
            frame1 = pktframe
            bat1 = pktbat
            hbat1 = hdrbat
            currint = None

        if pktframe < frame1:
            continue

        # TODO: pktframes incremetn by < 1 integration over FPGAs. Grrr.
        # Remove INT when we're able
        intno = ((pktframe - frame1)/intime)
        #print 'INTNO', intno, pktframe, frame1, intime, pktframe-frame1, type(pktframe), type(frame1)
        intno_bat = (pktbat - bat1) / batint

        if currint is None:
            currint = intno
            vints[currslot] = intno
            npkts = 0

        last_pktbat, last_pktframe, last_hdrbat, last_seqno, last_pktdate  = last_times.get((address, dev), (0, 0, 0, 0, None))
        datazero = np.all(data == 0)
        datazero_str = 'Z' if datazero else 'F'
        if last_pktdate is not None and pktdate is not None:
            pktdiff = str((pktdate - last_pktdate).microseconds/1e3)
            pktdatestr = str(pktdate.time())
        else:
            pktdiff = ''
            pktdatestr = ''

        
        if values.verbose:
        
            print('\t'.join(map(str, ('%s/%s' % (address[0], hdr['srcDevice']),
                                  pktdatestr, pktdiff,
                                  seqno, seqno-last_seqno,
                                  hdrbat, hdrbat - hbat1, hdrbat - last_hdrbat,
                                  pktbat, pktbat-bat1,  pktbat-last_pktbat, '%0.3f' % intno_bat,
                                  pktframe, pktframe - frame1,pktframe - last_pktframe, '%0.3f' % intno, start_frames[1:] - start_frames[0],
                                  'fdiff', fdiff, bdiff, chan, '%d%s' % (len(rawdata), datazero_str)
                                  ))))

        last_times[(address, dev)] = (pktbat, pktframe, hdrbat, seqno, pktdate)

if __name__ == '__main__':
    _main()
