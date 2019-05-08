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
import struct
from astropy.time import Time

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-o','--outfile', help='Output file', type=FileType('w'))
    parser.add_argument('-w','--max-boxcar', help='max width to declare', default=32, type=int)
    parser.add_argument('-d','--min-dm', help='minimum dm to declare', default=0, type=float)
    parser.add_argument('-b','--min-sn', help='minimum S/N to declare', default=0, type=float)
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
        # Ryan's addition
        sock.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR, 1)
        sock.bind((host, port))

        # add to multicast group # https://stackoverflow.com/questions/603852/how-do-you-udp-multicast-in-python
        mreq = struct.pack('4sl', socket.inet_aton(host), socket.INADDR_ANY)
        sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)
        total_cand = 0
        while True:
            data, addr = sock.recvfrom(1500)
            npdata = np.fromstring(data, sep=' ')
            npdata.shape = (-1, 8)
            total_cand += npdata.shape[0]
            sn = npdata[:, 0]
            sampnum = npdata[:, 1]
            tsec = npdata[:, 2]
            width = npdata[:, 3]
            idt = npdata[:, 4]
            dm = npdata[:, 5]
            beamno = npdata[:, 6]
            mask = (sn > values.min_sn) & (width < values.max_boxcar) & (dm > values.min_dm) & (beamno != 35) & (tsec > 10.)
            if np.any(mask):
                goodat = npdata[mask, :]
                good_sn = goodat[:, 0]
                best_cand_idx = np.argmax(good_sn)
                best_cand = goodat[best_cand_idx, :]
                sn, sampnum, tsec, width, idt, dm, beamno, cand_mjd = best_cand
                #print 'FOUND CANDIDATE', np.array2string(best_cand, precision=1), best_cand.shape
                now = Time.now()
                latency_ms = (now.mjd - cand_mjd)*86400.0*1e3
                best_beam = int(beamno)
                s = 'Found CANDIDATE: sn={} width={} dm={} mjd={} latency={}ms beam={} total_cands={}'.format(sn, width, dm, cand_mjd, latency_ms, best_beam, total_cand)
                print s
                cand_list = list(best_cand)
                cand_list.append(latency_ms)

                outf = values.outfile
                if outf:
                    outf.write('# S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, latency_ms\n')
                    outf.write('{:0.2f} {:0.0f} {:0.4f} {:0.0f} {:0.0f} {:0.2f} {:0.0f} {:0.10f} {:0.2f}\n'.format(*cand_list))
                    outf.write('# total_cand: {}\n'.format(total_cand))
                    outf.flush()
                    outf.close()
                    
                sys.exit(best_beam+10)
                
if __name__ == '__main__':
    _main()
