#!/usr/bin/env python
"""
CRAFT beamformer wrapper

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'


# import register definitions
from askap.adbe.regdefs import *

# import the Beamformer object from the library
from askap.adbe import Beamformer, adbe_msg_str, FilterTypes, ZoomModes, AdbeStatus, EventData1D, EventData, RBType

import atexit
import sys
import struct
from collections import namedtuple


NBEAMS = 72
NCHANS = 8




class NamedStruct(object):
    def __init__(self, name, fmt, attnames):
        self.struct = struct.Struct(fmt)
        self.named_tuple = namedtuple(name, attnames)
        
    def make_from_buffer(self, data):
        return self.named_tuple._make(self.struct.unpack(data[0:self.struct.size]))

    def pack(self):
        return self()

    def __call__(self, *args, **kwargs):
        t = self.named_tuple(*args, **kwargs)
        b = self.struct.pack(*t)
        return b
        

RxHeaderV1 = NamedStruct('RxHeaderV1', '<BBHBBBBH', 'ver packetSequence fragmentID packetType port packetID flags status')

RxHeaderV2 = NamedStruct('RxHeaderV2', '<BBHBBBBHBBI', 'ver packetSequence fragmentID packetType srcDevice packetID flags status srcPort spare bat')

TxHeader = NamedStruct('TxHeader', '<BBHBBHBHBI', 'ver packetSequence spare coreCommand port packetLength device wordCount command address')


def buffer_to_header(data):
    hdrv1 = RxHeaderV1.make_from_buffer(data)
    if hdrv1.ver == 0x2 and len(data) >= RxHeaderV2.struct.size:
        hdrtype = RxHeaderV2
        hdr = RxHeaderV2.make_from_buffer(data)
    else:
        hdrtype = RxHeaderV1
        hdr = hdrv1

    payload = data[hdrtype.struct.size:]
    return hdr, payload

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
        print(self.ctrl_addr, self.data_addr, rbtype)
        # NB: Beamformer.cc Boost wrapper doesn't have the port number in it
        Beamformer.__init__(self, self.ctrl_addr, self.data_addr, rbtype)

        atexit.register(self.check_stop_craft)

    def connect(self):
        code = Beamformer.connect(self)
        if code != AdbeStatus.ADBE_OK:
            raise ValueError('Couldnt connect')

        return code

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


    def wait_for_done(self):
        bf = self
        bf.getBullantCard().flushEvents()
        bf.getRedbackCard().clearPushFlags(0x3f, 0x20)
        timeout = 0
        while (timeout < 10):
            start = bf.getBullantCard().getBatFrame()
            if start.batStatus & 0x20 == 0: 
                break
            timeout += 1
                
            if timeout == 10:
                print('Timout - no BAT captured')
                return

    def start_craft(self):
        card = self
            
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

        card.getRedbackCard().clearPushFlags(0x3f, 0x20)

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
