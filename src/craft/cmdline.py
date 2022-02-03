#!/usr/bin/env python
"""
Utilities for parsing the command line

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'

import argparse

def strrange(rangestr):
    ''' Converts a range string into a list containing each value in that range
    the range strine can be comma separated integers, or hyphen separated integers 
    (endpoints included) or both.

    :rangestr: Input range as a string
    :see: argparse python module
    :raises: argparser.ArgumentTypeError if rangestr is badly formatted

    Examples:
    
    >>> strrange('1')
    [1]
    
    >>> strrange('1,2')
    [1, 2]
    
    >>> strrange('1-3')
    [1, 2, 3]
    
    >>> strrange('1-3,4,7,9-11')
    [1, 2, 3, 4, 7, 9, 10, 11]

    '''
    fullrange = []
    for c in rangestr.split(','):
        cbits = c.split('-')
        if len(cbits) == 1:
            fullrange.append(int(cbits[0]))
        elif len(cbits) == 2:
            start, end = list(map(int, cbits))
            fullrange.extend(list(range(start, end+1)))
        else:
            raise argparse.ArgumentTypeErrror("Invalid range string %s" % (rangestr))
    
    return fullrange


class Stime(object):
    
    def __init__(self, s):
        dt = None
        tint = None
        tfloat = None
        try:
            dt = dateutil.parser.parse(s)
        except ValueError: # Invalid date
            try:
                tint = int(s)
            except ValueError: # Not an int
                try:
                    tfloat = float(s)
                except ValueError: # not a float either
                    raise ValueError('Invalid time {}'.format(s))
    
