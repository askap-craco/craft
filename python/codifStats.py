#!/usr/bin/env python3

import sys, struct
import numpy as np


class CodifHeader:
    def __init__(self, bytes):
        (seconds, frame,framelength, antid, nchan, groupid, period,
         reserved1, numsamples, sync, reserved2) = struct.unpack('<IIIIIIIIQII', bytes[:48])

        self.invalid = seconds>>31
        self.complex = (seconds>>30)&0x1
        self.seconds = seconds & 0x3FFFFFFF

        self.nbits = ((framelength>>24)&0x1F)
        self.version = framelength >> 29
        self.framelength = (framelength & 0xFFFFFF)*8

        self.refepoch = (antid>>26)&0x3F
        self.representation = (antid>>22)&0xF
        self.antid = antid & 0xFFFF

        self.sampleblocklength = (nchan>>24)&0x1F
        self.nchan = nchan & 0xFFFF

        self.threadid = (groupid>>16)&0xFFFF
        self.groupid = groupid & 0xFFFF

        self.period = period & 0xFFFF

if len(sys.argv)<2:
    print("Usage: codifStats.py <file> [<file> ...]")


# Histogram 4 bit complex data
def frameHisto(frame, nchan):
    histo = np.zeros((nchan*2,16), dtype=int)
    nbyte = len(frame)
    assert (nbyte % nchan)==0, "Frame size must be multiple of # channels"

    for i in range(nbyte//nchan):
        for j in range(nchan):
            ireal = frame[i*nchan+j]&0xF;
            icomplex = frame[i*nchan+j]&0xF;
            histo[j*2][ireal] +=1;
            histo[j*2+1][icomplex] +=1;

    return(histo)


maxFrame = 10

for file in sys.argv[1:]:

    nframe = 0
    with open(file, "rb") as f:
        header = CodifHeader(f.read(64))

        iscomplex = header.complex
        bits = header.nbits
        nchan = header.nchan
        framelength = header.framelength

        histo = np.zeros((header.nchan*2,16), dtype=int)

        while nframe<maxFrame:
            data = f.read(framelength)
            if not data: break
            nframe += 1
            histo += frameHisto(data, nchan)

            h = f.read(64)
            if not h: break

    shuffleIdx = list(range(8,16))+list(range(8))
    for n in range(-8,8): print("{:5d} ".format(n),end="")
    print()
    
    for r in histo[:,shuffleIdx]:
        sum = r.sum()
        norm = r.astype(float)/sum*100
        for n in norm: print(" {:5.2f}".format(n),end="")
        print()
