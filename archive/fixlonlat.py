#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import os
import sys
import logging
import json

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    for fname in values.files:
        foutname = fname+'v2'
        fin = open(fname, 'rU')
        fout = open(foutname, 'w')
        for line in fin:
            s = json.loads(line)
            if 'index' in s.keys():
                fout.write(line)
            else:
                s['ant_direction'] = (s['ant_direction'][1], s['ant_direction'][0])
                s['beam_directions'] = [(ra, dec) for (dec, ra) in s['beam_directions']]
                for b in s['filterbanks']:
                    b['beam_direction'] = b['beam_direction'][1], b['beam_direction'][0]

                fout.write(json.dumps(s))
                fout.write('\n')
    

if __name__ == '__main__':
    _main()
