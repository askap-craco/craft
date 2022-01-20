#!/usr/bin/env python
"""
Emulates the CRAFT functionality in a beamformer

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys

# import register definitions
from askap.adbe.regdefs import *

# import the Beamformer object from the library
from askap.adbe import Beamformer, adbe_msg_str, FilterTypes, ZoomModes, AdbeStatus, EventData1D, EventData

from askap.craft.sigproc import SigprocFile
import askap.craft.pktdump as pktdump
import askap.time
import pickle as pickle
import atexit
import numpy as np
import time
import socket
import struct
import pylab
import pytz
import datetime
import threading
import select

from .craft_fil import TxHeader, RxHeaderV2


NBEAMS = 72
NCHANS = 8
NFPGAS = 6
SAMP_RATE = 32./27.*1e6 # samples/sec


class Emulator(object):
    def __init__(self, port,hdr=None, bfnum=None, pktfile_name=None):
        self.tx_addr = None
        self.is_streaming = False
        self.hdr = hdr
        self.socks = []
        self.tx_addr = {} # maps beamformer number to TX address
        self.bfnum = bfnum
        self.sock_map = {}
        
        if hdr is None:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            self.sock.bind(("127.0.0.1", port))
            assert self.bfnum is not None
            self.socks = [self.sock]
            self.int_cycles = 3
            self.int_time = 1000
        else:
            self.bw = float(hdr['BW'])
            self.int_time = int(hdr['INT_TIME'])
            self.int_cycles = int(hdr['INT_CYCLES'])

            if self.bfnum is None:
                bfs = np.arange(int(hdr['NUM_BEAMFORMERS']))
            else:
                bfs = [bfnum]
                
            self.freqs = []
            self.chans = []
            self.ports = []
            logging.debug('Makeing beamforer sockets for  beamforemrs: (%s), %s', len(bfs), bfs)
            for bfnum in bfs:
                addr = hdr['BEAMFORMER%d_ADDR' % bfnum]
                f = np.array(list(map(float, hdr['BEAMFORMER%d_FREQMAP' % bfnum].split(',')))) # MHz
                c = list(map(float, hdr['BEAMFORMER%d_CHANMAP'% bfnum].split(',')))
                p = int(hdr['BEAMFORMER%d_PORT' % bfnum])
                sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
                rxaddr = ("127.0.0.1", p)
                logging.debug('Binding to socket %s', rxaddr)
                sock.bind(rxaddr)
                self.freqs.append(f)
                self.chans.append(c)
                self.ports.append(p)
                self.socks.append(sock)
                self.sock_map[addr] = sock
                logging.info('bfnum %s with adresss %s listining on port %s', bfnum, addr, p)

            logging.debug('Got %s sockets', len(self.socks))

        if self.int_cycles > 7:
            raise NotImplemented('I cannnot be bothered doing fragmenting for large packets')

        self.pktfile_name = pktfile_name
        if self.pktfile_name is not None:
            self.pktfile = pktdump.pktopen(self.pktfile_name)

        self.sleep_time = float(self.int_time*self.int_cycles)/SAMP_RATE
        
    def process(self):
        logging.debug('ABout to select with %s sockets', len(self.socks))
        rsocks, wsocks, xsocks = select.select(self.socks, [], [])
        sock = rsocks[0]

        print('Got data from socket', sock, self.socks, rsocks)

        data, addr = sock.recvfrom(9000)
        hdr = TxHeader.make_from_buffer(data)

        logging.debug("got packet %s", str(hdr))
    
        if hdr.ver != 0x1:
            logging.info('Unknown packet version:{} '.format(hdr))
            return None

        if hdr.coreCommand == 0x80:
            # this means start streaming to the source address
            bfno = self.socks.index(sock) + 1
            assert bfno > 0
            self.tx_addr[bfno] = addr

            logging.info('Sending data for bf %s to address %s', bfno, addr)
            # send reply- no sure if this is correct, but craft_fil needs it

            if not self.is_streaming: 
                self.start_streaming() # sets up start_bat and other useful timestamps

            logging.info("Starting streaming %s %s %s %s", self.tx_addr, self.start_utc, self.start_bat, self.sleep_time)
            rxhd = RxHeaderV2(ver=0x2, packetSequence=0, fragmentID=0, packetType=128, srcDevice=0, packetID=0,
                              flags=0, status=0, srcPort=0, spare=0, bat=self.curr_bat&0xffffffff)

            logging.info('Sending reply %s bytes', len(rxhd))

            sock.sendto(rxhd, addr)


    def stop_streaming(self):
        logging.info("Stopping streaming")
        self.is_streaming = False
        self.streaming_thread.join()
        logging.info('Streaming stopped')

    def start_streaming(self):
        self.start_utc = datetime.datetime.now(tz=pytz.utc)
        self.start_bat = askap.time.utcDt2bat(self.start_utc)
        self.start_frameno = 0
        self.bat_increment = int(self.int_time*32./27.)
        self.frame_increment = self.int_time
        self.curr_frameno = self.start_frameno
        self.curr_bat = self.start_bat
        self.intnum = 0
        self.streaming_thread = threading.Thread(target=self.push_packets)
        self.streaming_thread.daemon = False
        self.is_streaming = True
        self.streaming_thread.start()

    def _make_times(self):
        ic = self.int_cycles
        times = np.zeros(4*ic, dtype=np.uint64)
        for cyc in range(self.int_cycles):
            times[2*cyc+0] = self.curr_frameno # start frame
            times[2*cyc+1] = self.curr_frameno + self.int_time - 1 # stop frame
            self.curr_frameno += self.frame_increment
            
            times[2*cyc+ic*2+0] = self.curr_bat # start bat
            times[2*cyc+ic*2+1] = self.curr_bat + self.bat_increment # stop bat
            self.curr_bat += self.bat_increment

            self.intnum += 1


        self._times = times
        #print 'FRAMEstart', self._times[0], self.frame_increment

        return times

    def get_data_noisonly(self, fpga):
        ic = self.int_cycles

        data = np.random.randn(NCHANS, ic, NBEAMS)+10.
        return data

    def get_data_pulses(self, fpga, pulse_interval=15000):
        ic = self.int_cycles

        data = np.random.randn(NCHANS, ic, NBEAMS)+128.
        inum = self.intnum
        minf = min(self.freqs)/1e3 # GHz
        maxf = max(self.freqs)/1e3 # GHz

        tsamp = self.int_time / SAMP_RATE # seconds
        tsamp_ms = tsamp*1e3

        for b in range(NBEAMS):
            dm = 10.*float(b) 

            for t in range(ic):
                tnow = inum + t # integrations
                ipulse = tnow / pulse_interval
                tpulse_start = ipulse*pulse_interval
                
                for c in range(NCHANS):
                    c2 = fpga*NCHANS + c
                    f = self.freqs[c2]/1e3 # GHz
            
                    tau = 4.15*dm*(f**-2 - maxf**-2) # milliseconds
                    assert tau >= 0, 'Invalid tau: %f'%tau
                    tau_int = tau / tsamp_ms
                    assert tsamp_ms >= 0
                    
                    x = tpulse_start + tau_int - tnow
                    inpulse = abs(x) < 1.0
                    
                    if inpulse:
                        data[c, t, b] += 5
#                        if b == 10:
#                            print tsamp, tsamp_ms, b, t, dm, tnow, tpulse_start, c,  c2, f, tau, tau_int, x, inpulse

        return data




    def get_data_ramps(self, fpga):
        ic = self.int_cycles
        im = self.intnum
        data = np.zeros((NCHANS, ic, NBEAMS))
        for c in range(NCHANS):
            for t in range(ic):
                tmod = (im + t) % 256
                data[c, t, :] = np.arange(NBEAMS) + tmod + self.freqs[fpga*NCHANS + c]

        return data
        

    def make_packet(self):
        ic = self.int_cycles
        times = self._make_times()
        for fpga in range(NFPGAS):
            # According to beamformer.cc it should be FTB with beam fastest
            data = self.get_data_pulses(fpga)
                #data = self.get_data_ramps(fpga)
            data = data.astype(np.float32).flatten()
            fpgaid = 1<<fpga
            rxhd = RxHeaderV2(ver=0x2, packetSequence=0, fragmentID=0, packetType=129, srcDevice=fpgaid, packetID=0,
                              flags=0, status=0, srcPort=0, spare=0, bat=self.curr_bat&0xffffffff)
            buf = rxhd + times.tostring() + data.tostring()


        # TODO: Work out which address to return
        return buf, None

    
    def make_replay_packet(self):
        assert self.pktfile is not None
        
        while True:
            hdr, address, data, dt  = next(self.pktfile)
            rxhd = RxHeaderV2(**hdr)
            cardno = int(address[0].split('.')[3])
            if cardno in list(self.tx_addr.keys()):
                break;


        allbytes = rxhd + data

        return allbytes, cardno

    def push_packets(self):
        ic = self.int_cycles
        try:
            while self.is_streaming:
                if self.pktfile_name is None:
                    buf, cardno = self.make_packet()
                else:
                    buf, cardno = self.make_replay_packet()

                sock = self.socks[cardno-1]
                addr = self.tx_addr[cardno]
                sock.sendto(buf, addr)
                time.sleep(self.sleep_time)
        except:
            logging.exception('Finished push pckets')


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Captures CRAFT autocorrelation spectra and writes to SIGPROC file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-i','--int-time', help='Integration time (1-65535)', default=1000, type=int)
    parser.add_argument('-s','--int-cycles', help='Number of cycles to combine (1-3)', default=3, type=int)
    parser.add_argument('-p','--port', help='UDP Listen port', default=12290, type=int)
    parser.add_argument('-H','--header', help='Header file to setup on')
    parser.add_argument('file', help='Replay file')
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    hdr = None
    if values.header is not None:
        hdr = {}
        for line in open(values.header, 'rU'):
            name, value = line.split(None, 1)
            hdr[name.strip()] = value.strip()

    int_time = int(values.int_time)
    assert int_time >=1 and int_time<= 65535, 'Invalid integration time {}'.format(int_time)
    
    int_cycles = int(values.int_cycles)
    assert int_cycles >=1 and int_cycles <= 7, 'Invalid int cycles {}'.format(int_cycles)

    emu = Emulator(values.port, hdr, None, values.file)
    while True:
        try :
            emu.process()
        except KeyboardInterrupt:
            sys.exit(0)


if __name__ == '__main__':
    _main()
