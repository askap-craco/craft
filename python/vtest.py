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
import dada
from crafthdr import DadaHeader
import datetime
import astropy
from astropy.time import Time
import aktime
import pytz

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

mode_to_nbits = (16, 8, 2, 1)
SAMP_RATE = 1e6*32./27.
NCHANS = 8

def coherent_disperse():
    h = np.exp(2*pi*1j*d/(f + f0))


def _main():
    from argparse import ArgumentParser
    import argparse
    parser = ArgumentParser(description='Script description')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-m','--mode', help='CRAFT Mode', type=int, default=0)
    parser.add_argument('-n','--nsamps', help='Number of samples', type=int, default=1024)
    parser.add_argument('-b','--beamid', help='Beam ID downloaded', type=int, default=0)
    parser.add_argument('-a','--antenna', help='Antenna number', type=int, default=1)

    parser.add_argument(dest='file', help='Output file', type=argparse.FileType('w'))
    parser.add_argument('--ascii-output', help='Output file as ascii', type=argparse.FileType('w'), default=None)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    hdr = make_header(values)
    print hdr
    d = np.random.randn(values.nsamps, NCHANS) + 1j*np.random.randn(values.nsamps, NCHANS)
    d *= 2*14
    # INJECT SIGNALS HERE ON A PER CHANNEL BASIS
    hdr.tofile(values.file)
    enc_func = mode_funcs[values.mode & 0x3]
    dout = enc_func(d)
    dout_bytes = dout.tostring()
    values.file.write(dout_bytes)
    values.file.close()
    if values.ascii_output:
        values.ascii_output.write(str(hdr))
        dri = np.zeros((values.nsamps, NCHANS, 2))
        dri[:,:,0] = np.real(d)
        dri[:,:,1] = np.imag(d)
        np.savetxt(values.ascii_output, d.flat)

def encode_16b16b(d):
    '''
    Encodes data as 16 bit signed integers

    >>> encode_16b16b(np.arange(-3, 3).reshape(6,1))
    
    '''
    nsamps, nchans = d.shape
    encdata = np.zeros((nsamps, nchans, 2), dtype=np.int16)
    encdata[:, :, 0] = np.round(np.real(d))
    encdata[:, :, 1] = np.round(np.imag(d))

    return encdata

def encode_8b8b(d):
    nsamps, nchans = d.shape
    encdata = np.zeros((nsamps, nchans, 2), dtype=np.int8)
    encdata[:, :, 0] = np.round(np.real(samps))
    encdata[:, :, 1] = np.round(np.imag(samps))

    return encdata

def encode_4b4b(d):
    '''
    Packs real & imag into 2 4-bit numbers
    '''
    nsamps, nchans = d.shape
    #encdata = np.zeros((nsamps, nchans), dtype=np.uint8)
    r = np.round(np.real(samps)).astype(uint8)
    i = np.round(np.imag(samps)).astype(uint8) & 0xf

    encdata = r + (i << 4)
    return encdata.tostring()

def encode_1b1b(d):
    nsamps, nchans = d.shape
    raise NotImplemented('Sorry, just ran out of time')
                       
mode_funcs = [encode_16b16b, encode_8b8b, encode_4b4b, encode_1b1b]


def make_header(values):
    mode = values.mode
    nbits_per_samp = mode_to_nbits[mode & 0x3] * 2
    nwords = values.nsamps / nbits_per_samp
    fpga_id = 1
    cardno = 2
    freqs = np.linspace(1200, 1208, NCHANS, endpoint=True)
    assert len(freqs) == NCHANS
    now = datetime.datetime.utcnow()
    now = now.replace(tzinfo=pytz.UTC)
    now_bat = aktime.utcDt2bat(now)
    now_time = Time(now, scale='utc')

    start_bat = now_bat - 2000000 # 2 seconds
    end_bat = start_bat + int(values.nsamps/SAMP_RATE*1e6)
    start_bat32 = now_bat  & 0xffffffff
    finish_bat32 = end_bat & 0xffffffff
    trigger_bat32 = start_bat32
    start_frameid = 1234
    finish_frameid = start_frameid + values.nsamps
    trigger_frameid = start_frameid

    h = DadaHeader()
    h += 'HDR_SIZE', 8192, 'Header size'
    h += 'DATA_TYPE', 'CRAFT_VOLTAGES', 'CRAFT_VOLTAGES'
    h += 'CMD_ARGS',  ' '.join(sys.argv), 'Command line arguments'
    h += 'SAMP_RATE', SAMP_RATE, 'Sample rate in samples per second'
    h += 'CRAFT_MODE', values.mode, 'The craft mode which describes the applicable bit depth (hex) among other things'
    h += 'NWORDS', nwords, 'The number of words downloaded'
    h += 'NSAMPS', values.nsamps, ' The number of samples downloaded'
    h += 'NBITS', nbits_per_samp, 'Number of bits per complex sample'
    h += 'OVERSAMPLING', '32/27', 'Oversampling factor (N/M)'
    h += 'BEAMID', values.beamid, 'The beam ID requested for the download *0-71, I think*'
    h += 'START_WRITE_BAT32', hex(start_bat32), 'Lowest 32 bits of BAT of the start of the buffer (hex) - direct from hardware'
    h += 'FINISH_WRITE_BAT32', hex(finish_bat32), 'Lowest 32 bits of BAT of the end of the buffer (hex)'
    h += 'TRIGGER_BAT32', hex(trigger_bat32), 'BAT   lowest 32 bits of BAT when the trigger occurred (hex)'
    h += 'START_WRITE_FRAMEID', start_frameid, 'Frame ID of the start of the buffer (decimal)'
    h += 'FINISH_WRITE_FRAMEID', finish_frameid, 'Frame ID of the end of the buffer (decimal)'
    h += 'TRIGGER_FRAMEID', trigger_frameid, 'Frame When the trigger occurred (decimal)'
    h += 'FPGA_ID', fpga_id, 'The FPGA from which this data was downloaded'
    h += 'CARD_NO', cardno, 'The card number from which this data was downloaded'
    h += 'ANT', 'AK%02d' % values.antenna, 'Antenna name (AKNN)'
    h += 'ANTENNA_NO', values.antenna, 'The antenna number from which this data was downloaded 1-36'
    h += 'NCHANS', NCHANS, 'Number of frequency channels'
    h += 'FREQS', ','.join(map(str,freqs)), 'Comma separated list of frequencies applying to each of channels in this file(MHz)'
    h += 'UTC_NOW', now.isoformat(), 'UTC date stamp for when the file was written (ISO format)'
    h += 'MJD_NOW', now_time.mjd, 'MJD date stamp for when the file was written'
    h += 'BAT_NOW', hex(now_bat), 'Current BAT all bits (hex)'

    h.reset_hdr_size()

    return h

    

if __name__ == '__main__':
    _main()
