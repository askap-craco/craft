#!/usr/bin/env python
"""
Listens to fredda on a udp socket

Copyright (C) CSIRO 2015
"""
import numpy as np
import os
import sys
import logging
import socket

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='hostport', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for p in values.hostport:
        host, port = p.split(':')
        port = int(port)
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        sock.bind((host, port))
        while True:
            data, addr = sock.recvfrom(1500)
            npdata = np.fromstring(data, sep=' ')
            npdata.shape = (-1, 7)
            sn = npdata[:, 0]
            sampnum = npdata[:, 1]
            tsec = npdata[:, 2]
            width = npdata[:, 3]
            idt = npdata[:, 4]
            dm = npdata[:, 5]
            beamno = npdata[:, 6]
            mask = (sn > 10) & (width < 10) & (dm > 30) & (beamno != 35)
            if np.any(mask):
                goodat = npdata[mask, :]
                print 'FOUND CANDIDATE'
                print goodat
                break



if __name__ == '__main__':
    _main()
