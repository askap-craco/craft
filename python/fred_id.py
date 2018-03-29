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
import sigproc as sgp
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
'''
dev:python3

python2 may have int/float bugs in this code, plz see comments below
'''
__author__ = "CRAFT Harry Qiu <hqiu0129@physics.usyd.edu.au>"

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    #parser.add_argument('-s','--show', action='store_true', help='Show')
    parser.add_argument('--probe', action='store_true', help='Show')
    parser.add_argument('--noname', action='store_true', help='Show')
    parser.add_argument('-o','--output',type=str,default='short_csv')
    parser.add_argument('-c','--column',type=str,default='name p0 dm raj decj',help='')
    parser.add_argument('-x','--threshold',type=float,default=10,help='fredda threshold')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    for i in values.files:
        filterbank=i
        ra,dec,dm_list=fredda_read(filterbank,fred_check=values.noname)
        print(ra,dec)
        print("Searching "+str(len(dm_list))+" Candidates")
        searcher(filterbank,ra,dec,dm_list)

    ##### search psrcat
def fredda_read(filterbank,fred_check=False):
    ##### read in filterbank sigproc format
    readin=sgp.SigprocFile(filterbank)
    fredda_file=filterbank+".cand.fof"
    if fred_check:
        fredda_file="fredda.cand.fof"
    ra0=readin.header['src_raj']
    dec0=readin.header['src_dej']
    if dec0 >= 0:
        pm=1
        pmsign=' +'
    else:
        pm=-1
        dec0=-dec0
    ####parse the ra and dec
    ra_hr="{0:02}".format(int(ra0//10000))
    ra_min="{0:02}".format(int((ra0-10000*(ra0//10000))//100))
    ra_sec="{0:02}".format((ra0-10000*(ra0//10000)-100*((ra0-10000*(ra0//10000))//100)))
    de_hr="{0:02}".format(int(dec0//10000))
    de_min="{0:02}".format(int((dec0-10000*(dec0//10000))//100))
    de_sec="{0:02}".format((dec0-10000*(dec0//10000)-100*((dec0-10000*(dec0//10000))//100)))
    coords=SkyCoord(ra_hr+" "+ra_min+" "+ra_sec+pmsign+de_hr+" "+de_min+" "+de_sec, unit=(u.hourangle,u.deg),frame='fk5')
    #open fredda
    fredda=np.loadtxt(fredda_file)
    dm_values=fredda.T[5][np.where(fredda.T[0]>values.threshold)]
    ra = coords.ra.deg
    dec = coords.dec.deg
    return ra,dec,dm_values
def searcher(filterbank,ra,dec,dm_values):
    f=open(filterbank[:-4]+".psrdb",'w')
    for dm in dm_values:
        logic=" RAJ > "+str(ra-1)+" && "+" RAJ < "+str(ra+1)+ "&& DECJ > "+str(dec-1)+" && "+" DECJ < "+str(dec+1) +" && DM > "+str(dm-5)+" && DM < "+str(dm+5)
        os.system("psrcat -c 'name dm raj decj' -l '"+str(logic)+"' -o short_csv > fredda_temp_psrcat")
        m=open("fredda_temp_psrcat",'r')
        mtext=m.readlines()
        m.close()
        detect=len(mtext)-2
        if detect > 0:
            print("Pulsar Identified\n")
            print(mtext[2:])
            f.write(filterbank+"\n")
            for i in xrange(detect):
                f.write(mtext[i+2])
    f.close()

if __name__ == '__main__':
    _main()
