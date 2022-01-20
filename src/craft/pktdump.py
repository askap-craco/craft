#!/usr/bin/env python
"""
Packet dump serialisation classes

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'

import struct
from collections import namedtuple
from . crafthdr import DadaHeader
from . beamformer import RxHeaderV1, RxHeaderV2, TxHeader, NamedStruct, buffer_to_header
import socket
import datetime
import pickle as pickle

CraftSpectrum_struct = struct.Struct('<QQQQ576f')

SpecHeader = NamedStruct('SpecHeader', '<QQQQ', 'startFrame stopFrame startBAT stopBAT')


class PicklePacketFile(object):
    def __init__(self, filename):
        self.fin = open(filename, 'r')
        self.hdr = pickle.load(self.fin)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            pkt = pickle.load(self.fin)
            # usually RxHeaderV2, address (IP string , port), pktdata (bytes), datetime
            return pkt
        except EOFError:
            raise StopIteration

class RawPacketFile(object):
    ''' Written by craft_udpdb in raw format'''
    hdr_struct = struct.Struct('!III4sHI')
    def __init__(self, filename, addr_swap=True):
        fin = open(filename, 'r')
        dada_hdr = DadaHeader.fromfile(filename)
        self.hdr = {}
        self.dada_hdr = dada_hdr
        for k, (v, comment) in dada_hdr.items():
            self.hdr[k] =v

        # Fix up some headers
        self.hdr['CRAFT_INT_CYCLES'] = int(self.hdr['INT_CYCLES'])
        self.hdr['CRAFT_INT_TIME'] = int(self.hdr['INT_TIME'])
        self.hdr['nbits'] = int(self.hdr['NBIT'])

        nbf = int(self.hdr['NUM_BEAMFORMERS'])
        cards = []
        for i in range(nbf):
            addr = self.hdr['BEAMFORMER{}_ADDR'.format(i)]
            card = addr.split('.')[3]
            cards.append(card)

        self.hdr['CARDS'] = ','.join(cards)

        hdr_nbytes = int(self.hdr['HDR_SIZE'])
        fin.seek(hdr_nbytes)
        
        self.addr_swap = addr_swap
        self.fin = fin

    def __iter__(self):
        return self

    def __next__(self):
        import binascii

        addr = '\x00\x00\x00\x00'
        while addr == '\x00\x00\x00\x00':
            hbuf = self.fin.read(RawPacketFile.hdr_struct.size)
            if len(hbuf) < RawPacketFile.hdr_struct.size:
                raise StopIteration()

            bits = RawPacketFile.hdr_struct.unpack(hbuf)

            (marker, sec, usec, addr, port, nbytes)  = bits
            
            if self.addr_swap:
                addr = addr[::-1]
                port = socket.ntohs(port)
        
            # hdr, address, data
            address = (socket.inet_ntoa(addr), port)

            dt = datetime.datetime.utcfromtimestamp(sec)
            dt = dt.replace(microsecond=usec)
            allbytes = self.fin.read(nbytes)
            if len(allbytes) < nbytes:
                raise StopIteration()

        assert address[0] != '0.0.0.0'
        hdr, payload = buffer_to_header(allbytes)

        return (hdr._asdict(), address, payload, dt)


def pktopen(filename):

    try:
        pktfile = PicklePacketFile(filename)
    except pickle.PickleError as EOFError:
        pktfile = RawPacketFile(filename)

    return pktfile
        
