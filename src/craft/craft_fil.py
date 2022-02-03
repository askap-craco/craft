#!/usr/bin/env python
"""
Downloads CRAFT data from a beamformer card

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys

# import register definitions
from askap.adbe.regdefs import *

# import the Beamformer object from the library
from askap.adbe import Beamformer, adbe_msg_str, FilterTypes, ZoomModes, AdbeStatus, EventData1D, EventData, RBType

from askap.craft.sigproc import SigprocFile
import askap.time
import pickle as pickle
import atexit
import numpy as np
import time
import socket
import struct
import pylab
import binascii
from collections import namedtuple

NBEAMS = 72
NCHANS = 8

class CraftBeamformer(Beamformer):
    '''
    Implents methods to handle events from 
    C++ library Beamformer object

    Events may come from multiple threads but Python wrapper
    ensures the calls to the handleXXX methods are serialized
    
    @see py-adbe/examples/Beamformer.py
    '''

    def __init__(self, ant, card):

        assert card >= 1 and card < 255
        assert ant >= 1 and ant < 255

        self.ctrl_addr = "10.1.%d.%d" % (ant, card)
        self.data_addr = "10.2.%d.%d" % (ant, card)
        self.console = CraftConsole()
        self.lastProgressUpdate = 0
        self._craft_running = False

        # for emulator, we increment port number
        # on same host so set incrementPortNum to True
        rbtype = RBType.RBTYPE_BEAMFORMER
        Beamformer.__init__(self, self.ctrl_addr, self.data_addr, rbtype)

        atexit.register(self.check_stop_craft)

    def callback(self,  msg):
        '''Required by adbe.Beamformer
        '''
        self.console.write("%s\n" % adbe_msg_str(msg.type))

    def showProgressUpdate(self, progress):
        '''Required by adbe.Beamformer

        display progress update
        '''
        if progress == 100:
            self.lastProgressUpdate = 0
            self.console.write('\n')
        else:
            if (progress > self.lastProgressUpdate):
                self.console.write('\b\b\b\b\b')
                self.console.write('x[%2.2d%%]' % progress)
                self.lastProgressUpdate = progress

    def handleInterrupt(self, card, value):
        ''' handle an interrupt from C++ library

        Required by adbe.Beamformer
        '''

        self.console.write('interrupt received card={} value={} \n'.format(card, value))


    def start_craft(self):
        card = self
        # Get ready for download
        card.getBullantCard().flushEvents()
        card.getRedbackCard().clearPushFlags(0x3f, 0x20)
            
        # Issue start
        eventsArray = EventData1D()
        event = EventData()
    
        event.eventTime = 0
        event.eventData = 0
        event.eventMask = 0xffff & ~(BF_EVENT_CRAFT_START)
        eventsArray.append(event)
        event.eventTime = 0
        event.eventData = BF_EVENT_CRAFT_START
        event.eventMask = 0xffff & ~(BF_EVENT_CRAFT_START)
        eventsArray.append(event)
        event.eventTime = 0
        event.eventData = 0
        event.eventMask = 0xffff & ~(BF_EVENT_CRAFT_START)
        eventsArray.append(event)
        card.getBullantCard().loadEvents(eventsArray)
        self._craft_running = True

    def stop_craft(self):
        card = self
        eventsArray = EventData1D()
        event = EventData()
        event.eventTime = 0
        event.eventData = 0
        event.eventMask = 0xffff & ~(BF_EVENT_CRAFT_STOP)
        eventsArray.append(event)
        event.eventTime = 0
        event.eventData = BF_EVENT_CRAFT_STOP
        event.eventMask = 0xffff & ~(BF_EVENT_CRAFT_STOP)
        eventsArray.append(event)
        event.eventTime = 0
        event.eventData = 0
        event.eventMask = 0xffff & ~(BF_EVENT_CRAFT_STOP)
        eventsArray.append(event)
        card.getBullantCard().loadEvents(eventsArray)

    def get_capability_register(self, fpga):
        assert 1 <= fpga <= 6

        d = self.getBullantCard().readRegister(fpga, DEVICE_REGS, int('103', 16))
        return d & 0xffffffff

    def check_craft_enabled(self):
        all_enabled = True
        for fpga in range(6):
            bits = self.get_capability_register(fpga+1)
            craft_enable = (bits & 0x10) == 0x10
            print('FPGA %d capabilities 0x%x CRAFT enabled? %s' % (fpga, bits, craft_enable))
            all_enabled &= craft_enable

        assert all_enabled, 'CRAFT not enabled on this beamformer'

        

    def check_stop_craft(self):
        if self._craft_running:
            self.stop_craft()

    def __del__(self):
        self.check_stop_craft()


class CraftConsole(object):
    def __init__(self):
        pass

    def run(self):
        pass

    def input(self, msg):
        return input(msg)

    def output(self, msg):
        sys.stdout.write(msg)
        sys.stdout.flush()

    write = output

class NamedStruct(object):
    def __init__(self, name, fmt, attnames):
        self.struct = struct.Struct(fmt)
        self.named_tuple = namedtuple(name, attnames)
        
    def make_from_buffer(self, data):
        return self.named_tuple._make(self.struct.unpack(data[0:self.struct.size]))

    def __call__(self, *args, **kwargs):
        t = self.named_tuple(*args, **kwargs)
        b = self.struct.pack(*t)
        return b
        

RxHeaderV1 = NamedStruct('RxHeaderV1', '<BBHBBBBH', 'ver packetSequence fragmentID packetType port packetID flags status')

RxHeaderV2 = NamedStruct('RxHeaderV2', '<BBHBBBBHBBI', 'ver packetSequence fragmentID packetType srcDevice packetID flags status srcPort spare bat')

TxHeader = NamedStruct('TxHeader', '<BBHBBHBHBI', 'ver packetSequence spare coreCommand port packetLength device wordCount command address')

CraftSpectrum_struct = struct.Struct('<QQQQ576f')
assert CraftSpectrum_struct.size == NBEAMS*NCHANS*4 + 32 #1184

SpecHeader = NamedStruct('SpecHeader', '<QQQQ', 'startFrame stopFrame startBAT stopBAT')

class CraftSpectrum(object):
    def __init__(self, b):
        s = SpecHeader.struct.size
        self.hdr = SpecHeader.make_from_buffer(b[0:s])
        self.data = np.fromstring(b[s:], dtype=np.float32)
        if len(self.data) == 72*8:
            self.data.shape = (72, 8)

assert TxHeader.struct.size == 16
assert RxHeaderV1.struct.size == 10
assert RxHeaderV2.struct.size == 16


class CraftPushReceiver(object):
    def __init__(self, hostport):
        self._hostport = hostport
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            
    def send(self, data, flags=0):
        self.sock.sendto(data, flags, self._hostport)

    def recv(self):
        data, address = self.sock.recvfrom(9000)
        hdrv1 = RxHeaderV1.make_from_buffer(data)
        if hdrv1.ver == 0x2 and len(data) >= RxHeaderV2.struct.size:
            hdrtype = RxHeaderV2
            hdr = RxHeaderV2.make_from_buffer(data)
        else:
            hdrtype = RxHeaderV1
            hdr = hdrv1

        return hdr, address, data[hdrtype.struct.size:]

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
    parser = ArgumentParser(description='Captures CRAFT autocorrelation spectra and writes to SIGPROC file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-a', '--ant', help='Antenna number', type=int, default=-1)
    parser.add_argument('-c','--card', help='Card number', type=int, default=-1)
    parser.add_argument('-o','--output', help='Output file name')
    parser.add_argument('-i','--int-time', help='Integration time (1-65535)', default=1000, type=int)
    parser.add_argument('-s','--int-cycles', help='Number of cycles to combine (1-7)', default=7, type=int)
    parser.add_argument('-b','--beam-number', help='Beam number to save voltages for (0-35). Unspecified for all beams', default=None)
    parser.add_argument('-m','--buffer-mode', help='Buffer number of bits', choices=[16, 8, 4, 1], type=int, default=16)
    parser.add_argument('--push-address', help='Push receiver address')
    parser.add_argument('--no-setup', dest='do_setup', action='store_false', help='Dont do craft setup')
    
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
    
    antno = int(values.ant)
    cardno = int(values.card)

    mjdnow = 55000.

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
           'ANTNO':antno,
           'CARDNO':cardno
           }



    if values.do_setup:
        bf = CraftBeamformer(antno, cardno)
        # Program beamformer setuff (nicely like)
        programFPGA = 0  # please dont' program FPGA - just be nice
        bandNo = 1 # 0 = 1200 MHz BPF, 1 = 1450 MHz BPF, 2 - 1800 MHz BPF, 3 - LPF
        centerFreq = 1140 # Centre frequency
        zoomMode = 0 # ZOOM_NONE
    
        if bf.connect() != AdbeStatus.ADBE_OK:
            raise ValueError("Couldn't connect")

        bf.check_craft_enabled()
        pushDownloadDelay = 0.0
        bf.startupBeamformer(programFPGA, FilterTypes.values[bandNo], int(centerFreq), ZoomModes.values[zoomMode], pushDownloadDelay)
        bf.getBullantCard().enableOffboardLinks()
    else:
        bf = None

    if values.push_address is None:
        hostport = ('10.2.%d.%d' % (antno, cardno), 12290)
    else:
        hostport = values.push_address.split(':')

    hostport = (hostport[0], int(hostport[1]))
    assert len(hostport) == 2
    print('Using push address', hostport)
    push_rx = CraftPushReceiver(hostport)
    print('Doing my own reset')
    push_rx.send_reset()

    print('And the result was')
    print(push_rx.recv())

    if values.do_setup:

        # setup CRAFT with specified stuff
        bf.setupCraft(int_time, int_cycles, bufmode, beam_number, 0)

        # and start that puppy up like a monkey
        bf.start_craft()
        bf.getRedbackCard().clearPushFlags(0x3f, 0x20)

    time.sleep(1)

    print('Recording...')

    pkl_file = open(values.output, 'w')
    pickle.dump(hdr, pkl_file)
    
    pylab.ion()
    #plots = Plots()
    while True:
        try:
            #result = bf.processCraftSpectra(2)
            hdr, address, data = push_rx.recv()
            if hdr.packetType != 129:
                continue

            print(hdr, address, len(data))#, binascii.hexlify(data)
            nint = int_cycles # get configured value, don't assume the packet size is right
            hdrsize = nint*4*8
            expected_size = hdrsize + 4*NBEAMS*NCHANS*nint
            actual_size = len(data)
            num_padbytes = expected_size - actual_size
            if expected_size != actual_size:
                print('Expected size', expected_size, 'actual size', actual_size)
            assert num_padbytes >= 0
            if num_padbytes > 0:
                data += '\x00'*num_padbytes

            times = np.fromstring(data[0:hdrsize], dtype=np.uint64)
            lt = len(times)
            assert lt % 2 == 0
            frames = times[0:lt/2]
            bats = times[lt/2:]
            startFrames = frames[0::2]
            stopFrames = frames[1::2]
            startBats = bats[0::2]
            stopBats = bats[1::2]


            print(hdr, 'BAT=', hex(hdr.bat), len(data), binascii.hexlify(data))
            print('\t startBATs', ' '.join(map(hex, startBats)))
            print('\t stopBATs', ' '.join(map(hex, stopBats)))
            print('\t diffs', (stopBats - startBats))
            print('\t startFrames', startFrames)
            print('\t stopFrames', stopFrames)
            print('\t diffs', (stopFrames - startFrames))

            spectra = np.fromstring(data[hdrsize:], dtype=np.float32)
            # According to beamformer.cc the it's in FTB ordering, with
            # beam going fastest
            spectra.shape = (NCHANS, nint, NBEAMS)
            
            #plots.update(spectra)
            
#            n = (1184*2/8)*8
#            du = np.fromstring(data[0:n], dtype=np.uint64)
#            for i in xrange(len(du)):
#                print i, hex(du[i])


        #if spfile is None:
        #spfile = SigprocFile(values.output, 'w', hdr)

            #write_result_pkl(pkl_file, result)
        except:
            logging.exception('Exception processing spectra')
            break
            
    pkl_file.close()
    if bf is not None:
        del bf

def write_result_pkl(pkl_file, results):
    nbeams = len(results)
    nints = len(results[0].point)
    nfreqs = len(results[0].point[0])
    data = np.empty((nbeams, nfreqs), dtype=np.float64)
    freqs = np.empty_like(data)
    fAs = np.empty((nbeams, nfreqs), dtype=np.uint64)
    fBs = np.empty_like(fAs)
    bAs = np.empty_like(fAs)
    bBs = np.empty_like(fAs)
    
    for i in range(nints):
        for beam in range(72):
            r = results[beam].point[i]
            data[beam,:] = [p.data for p in r]
            freqs[beam,:] = [p.freq for p in r]
            fAs[beam,:] = [p.startFrame for p in r]
            fBs[beam,:] = [p.stopFrame for p in r]
            bAs[beam,:] = [p.startBAT for p in r]
            bBs[beam,:] = [p.stopBAT for p in r]

        pickle.dump((i, data, freqs, fAs, fBs, bAs, bBs), pkl_file)
        
def bat2mjd(bat):
    utcdt = askap.time.bat2utcDt(bat)
    mjd = askap.time.utcDt2mjd(utcdt)

if __name__ == '__main__':
    _main()
