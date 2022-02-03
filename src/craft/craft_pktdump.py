#!/usr/bin/env python
"""
Downloads CRAFT data from a beamformer card

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys


from askap.craft.cmdline import strrange
from askap.craft.sigproc import SigprocFile
from askap.craft.beamformer import CraftBeamformer, TxHeader, buffer_to_header
import askap.time
import pickle as pickle
import atexit
import numpy as np
import time
import socket
import struct
import pylab
import binascii
import select
from collections import namedtuple
import datetime

from askap.adbe import adbe_msg_str, FilterTypes, ZoomModes, AdbeStatus, EventData1D, EventData, RBType

class CraftPushReceiver(object):
    def __init__(self, hostport):
        self._hostport = hostport
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            
    def send(self, data, flags=0):
        self.sock.sendto(data, flags, self._hostport)

    def recv(self):
        data, address = self.sock.recvfrom(9000)
        hdr, payload = buffer_to_header(data)
        return hdr, address, payload

    def send_reset(self):
        reset_hdr = TxHeader(ver=0x1, packetSequence=1, spare=0, coreCommand=0x80, port=0, packetLength=2, device=0, wordCount=0, command=0, address=0)
        self.send(reset_hdr)

import matplotlib.pyplot as plt
pylab.ion()
class Plots(object):

    def __init__(self):
        self.fig, self.axes2d = plt.subplots(3,2 , squeeze=True)
        self.axes = []
        for a in self.axes2d:
            self.axes.extend(a)

    def update(self, data):
        nint = data.shape[2]
        for i in range(min(len(self.axes), nint)):
            self.axes[i].imshow(data[i, :,:].T, aspect='auto')

        pylab.draw()


def _main():
    from argparse import ArgumentParser
    from askap.craft.cmdline import strrange
    parser = ArgumentParser(description='Captures CRAFT autocorrelation spectra and writes to SIGPROC file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-a', '--antennas', help='Antenna number, comma sparated or dash. e.g. "3", "13-15", "1,2,14"', type=strrange)
    parser.add_argument('-c','--cards', help='Card number. comma separated, or dash. e.g. "3", "1-7" or "1,2,4-7"', default="1-7", type=strrange)
    parser.add_argument('-o','--output', help='Output file name')
    parser.add_argument('-i','--int-time', help='Integration time (1-65535)', default=1000, type=int)
    parser.add_argument('-s','--int-cycles', help='Number of cycles to combine (1-7)', default=7, type=int)
    parser.add_argument('-b','--beam-number', help='Beam number to save voltages for (0-35). Unspecified for all beams', default=None)
    parser.add_argument('-m','--buffer-mode', help='Buffer number of bits', choices=[16, 8, 4, 1], type=int, default=16)
    parser.add_argument('--push-address', help='Push receiver address')
    
    parser.add_argument('--do-startup', dest='do_startup', action='store_true', help='Do beamformer startup')
    parser.add_argument('--do-setupcraft', dest='do_setupcraft', action='store_true', help='Do craft setup')

    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    int_time = int(values.int_time)
    assert int_time >=1 and int_time<= 65535, 'Invalid integration time {}'.format(int_time)
    
    int_cycles = int(values.int_cycles)
    assert int_cycles >=1 and int_cycles <= 7, 'Invalid int cycles {}'.format(int_cycles)

    samp_rate = 32./27.*1e6 # samples/sec
    tsamp = float(int_time)/samp_rate

    bufmode_map = {16:0, 8:1, 4:2, 1:3}
    assert values.buffer_mode in list(bufmode_map.keys()), 'Invalid buffer mode'
    bufmode = bufmode_map[values.buffer_mode]

    if values.beam_number is None:
        beam_number = -1
    else:
        beam_number = int(values.beam_number)
        assert beam_number >= 0 and beam_number < 36, 'Invalid beam number to save {}'.format(beam_number)
        bufmode += 4 # This tells the firmware to save only 1 beam
    
    ants = values.antennas
    cards = values.cards
    print('Cards', values.cards, cards)

    mjdnow = 55000.
    programFPGA = 0  # please dont' program FPGA - just be nice
    bandNo = 1 # 0 = 1200 MHz BPF, 1 = 1450 MHz BPF, 2 - 1800 MHz BPF, 3 - LPF
    centerFreq = 939. # Centre frequency
    zoomMode = 0 # ZOOM_NONE


    hdr = {'telescope_id':99, # Junk
           'TELESCOP':'ASKAP',
           'machine_id':99, # Junk
           'data_type':1, # filterbank
           'rawdatafile': values.output,
           'sourcename': 'test',
           'barycentric': 1,
           'pulsarcentric': 0,
           'az_start': 0., # azimuth at start of scan (degrees)
           'za_start': 0., # zenith at start of scan (degrees)
           'src_raj': 0.0, # RA (J2000) of source (hhmmss.s)??
           'src_dej': 0.0, # Dec (J2000) of source (ddmmss.s)??
           'tstart': mjdnow, # MJD of first sample
           'tsamp': tsamp, # Time (s) between samples
           'nbits': 32, # Number of bits per time sample
           'nsamples': 0, # number of time sampels in data file - dont use. File size is used instead
           'fch1': 710.0, # Centre frequency of first filterbank channel
           'foff': 1.0, # filterbank channel bandwidth
           'nchans': 48,
           'nifs': 72,
           'CRAFT_INT_TIME':int_time,
           'CRAFT_INT_CYCLES':int_cycles,
           'CRAFT_BUFMODE':bufmode,
           'CRAFT_BEAM_NO':beam_number,
           'ANTS':','.join(map(str, ants)),
           'CARDS': ','.join(map(str, cards)),
           'programFPGA':programFPGA,
           'bandNo':bandNo,
           'centerFreq':centerFreq,
           'zoomMode':zoomMode,
           'do_startup':values.do_startup
           }

# Program beamformer setuff (nicely like)
    bfs = []
    for antno in ants:
        bfs.extend([CraftBeamformer(antno, cardno) for cardno in cards])
    
    if values.do_startup:

        for bf in bfs:
            if bf.connect() != AdbeStatus.ADBE_OK:
                raise ValueError("Couldn't connect")

            bf.check_craft_enabled()
            pushDownloadDelay = 0.0
            bf.startupBeamformer(programFPGA, FilterTypes.values[bandNo], int(centerFreq), ZoomModes.values[zoomMode], pushDownloadDelay)
            bf.getBullantCard().enableOffboardLinks()


    if values.do_setupcraft:
        
        for bf in bfs:
            if bf.connect() != AdbeStatus.ADBE_OK:
                raise ValueError('Couldnt connect')

            # setup CRAFT with specified stuff
            bf.setupCraft(int_time, int_cycles, bufmode, beam_number, 0)
            bf.wait_for_done()

            # and start that puppy up like a monkey
            bf.start_craft()

    push_receivers = []
    sockets = []
    for antno in ants:
        push_addresses = [('10.2.%d.%d' % (antno, cardno), 12290) for cardno in cards]
        rxs = [CraftPushReceiver(addr) for addr in push_addresses]
        socks = [rx.sock for rx in rxs]

        push_receivers.extend(rxs)
        sockets.extend(socks) 
    
    for push_rx in push_receivers:
        print('Doing my own reset')
        push_rx.send_reset()
        
        print('And the result was')
        print(push_rx.recv())

    #time.sleep(1)

    print('Recording...')

    pkl_file = open(values.output, 'w')
    pickle.dump(hdr, pkl_file)

    assert len(sockets) == len(push_receivers)
    
    while True:
        try:
            rx, junk, junk = select.select(sockets, [], [])
            for rx_sock in rx:
                push_rx = push_receivers[sockets.index(rx_sock)]
                hdr, address, data = push_rx.recv()
                now = datetime.datetime.utcnow()
                pickle.dump((hdr._asdict(), address, data, now), pkl_file)
        except:
            logging.exception("Exception while receiving!")
            sys.exit(1)
               

if __name__ == '__main__':
    _main()
