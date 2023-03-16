#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
from craft.craco_wcs import *
from astropy.coordinates import Angle
from astropy import units as u


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def test_wcs2():
    target = SkyCoord('08h23m00s', '-46d00m00s')
    time = Time('2023-03-16T05:00:00')
    lmpix = Angle(['0.1d','0.1d'])
    Npix = 256
    inttime = 1.7*u.millisecond
    mywcs = CracoWCS(target, time, lmpix, Npix,inttime)

    assert mywcs.time == time
    assert mywcs.target == target
    print(mywcs.wcs2)
    print(mywcs.wcs2.to_header())
    print(mywcs.wcs3)
    #print(mywcs.wcs3.to_header())    # this is broken
    print(mywcs.wcs3_hdr)

    
    

if __name__ == '__main__':
    _main()
