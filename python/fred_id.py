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
#import sigproc as sgp
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
'''
dev:python3

python2 may have int/float bugs in this code, plz see comments below
'''
__author__ = "CRAFT Harry Qiu <hqiu0129@physics.usyd.edu.au>"


#### define psr/frb i/o operations
class FileOperations:
    def open_files(self,psr,frb):
        self.file1 = open(psr+'.psr', 'w')
        self.file2 = open(frb+'.frb', 'w')

    def write_psr(self,line):
        self.file1.write(line)
    def write_frb(self,line):
        self.file2.write(line)
    def close_files(self):
        self.file1.close()
        self.file2.close()
#remember to set class into function variables
#def func(a1,idio):
#    idio.write_psr(a1)
    ##### search psrcat
def read_fredda_hdr(file,x=39,y=40):
    infile=open(file,'r')
    hdr=infile.readlines()
    infile.close()
    sbid=hdr[9].split()[1]
    obsid=hdr[8].split()[1]
    ant=hdr[22].split()[1]
    ra=np.fromstring(hdr[x].split()[1],sep=',')[0:36]
    dec=np.fromstring(hdr[y].split()[1],sep=',')[0:36]
    pos=np.array([ra,dec]).T
    ######/data/TETHYS_1/craftop/auto_cands//SB01231_20180329103213_co18_tethys4.cand.fof.mbeam.2.good
    ####this isn't finished need fof file name format generation "SB{0:05}_".format(int(sbid))+" _co18_"+".cand.fof.mbeam.2.good"
    return pos
def fredda_read(fredda,beamno):
    '''
    ##### read in filterbank sigproc format
    readin=sgp.SigprocFile(filterbank)
    fredda_file=filterbank+".cand.fof"
    if fred_check:
        fredda_file='fredda.cand.fof'
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
    '''
    #open fredda
    #dm_values=fredda.T[5+1][np.where(fredda.T[6+1]==beamno)] ###+1 for good mode

    fof_array=np.where(fredda.T[7].astype(int)==beamno)[0]
    bm_values=fredda[fof_array]
    #ra = coords.ra.deg
    #dec = coords.dec.deg
    return bm_values,fof_array

def searcher(beamno,ra,dec,blist,foflines,idio,psrname,psrdm,psrra,psrdec):
    #searcher(i,ra,dec,dm_list,f1,f2,cat_name,cat_dm,cat_ra,cat_dec)
    #os.system("psrcat -c '"+format+"' > psrcat.csv")
    dm_values=blist.T[5+1] ##### dm values of candidates
    #select_ra=np.intersect1d(np.where(psrra>ra-1),np.where(psrra<ra+1))
    #select_dec=np.intersect1d(np.where(psrdec>dec-1),np.where(psrdec<dec+
    ##### query from psrcat array
    ### select pulsars within position
    select_pos=np.intersect1d(np.intersect1d(np.where(psrra>ra-5),np.where(psrra<ra+5)),np.intersect1d(np.where(psrdec>dec-5),np.where(psrdec<dec+5)))
    ##### extract DM information of selected pulsars
    posxdm=psrdm[select_pos]

    all_pos= np.where(psrra)[0]
    clen=len(posxdm)
    print(clen,'pulsars within beam')
    if clen>0:
        for i,match in enumerate(dm_values):
            print 'number',i
            writeline=foflines[i]
            print('dm',match)
            select_dm=np.intersect1d(np.where(posxdm>match-5),np.where(posxdm<match+5))
            if len(select_dm)>0:
                print 'Found Candidate'
                for j in select_dm:
                    real_pos=select_pos[j]
                    writename=psrname[real_pos]
                    print(writename,psrdm[real_pos])
                    idio.write_psr(writename+' '+writeline)
            else:
                idio.write_frb(writeline)

    else:
        for i in foflines:
            idio.write_frb(i)


    ####query from psrcat
    '''
    for i,dm in enumerate(dm_values):
        logic="RAJ > "+str(ra-1)+" && "+" RAJ < "+str(ra+1)+ "&& DECJ > "+str(dec-1)+" && "+" DECJ < "+str(dec+1) +" && DM > "+str(dm-5)+" && DM < "+str(dm+5)
        os.system("psrcat -c ' "+format+" ' -l '"+str(logic)+"' -o short_csv > fredda_temp_psrcat"+str(i))
    os.system("cat fredda_temp_psrcat* > psrcat_query")
    os.system("rm fredda_temp_psrcat*")
    '''

#def _main():
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
#parser.add_argument('-s','--show', action='store_true', help='Show')
parser.add_argument('--probe', action='store_true', help='Show')
parser.add_argument('--noname', action='store_true', help='Show')
#parser.add_argument('-o','--output',type=str,default='short_csv')
parser.add_argument('-f','--fof',type=str,default='') #### fredda candidate file input here
parser.add_argument('-c','--column',type=str,default='name p0 dm raj decj',help='')
parser.add_argument('-x','--threshold',type=float,default=10,help='fredda threshold')
parser.add_argument(dest='files', nargs='+') ####hdr file name here
parser.set_defaults(verbose=False)
values = parser.parse_args()
print values.files
if values.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
############ Generating psrcat database
if os.path.exists("psrcat.csv"):
    print("reading data from psrcat")
else:
    os.system("psrcat -c 'name p0 dm raj decj' -o short_csv -nohead -nonumber > psrcat_.csv ")
    print("generating psrcat.csv file")
    cat=np.genfromtxt("psrcat_.csv",skip_header=2,delimiter=';',dtype=str)
    name=cat.T[0]
    p0=cat.T[1]
    dm=cat.T[2]
    ra=cat.T[3]
    dec=cat.T[4]
    newcat=open("psrcat.csv",'w')
    for i in range(len(ra)):
         c = SkyCoord(ra[i]+' '+dec[i], unit=(u.hourangle, u.deg))
         newcat.write(name[i]+";"+p0[i]+";"+dm[i]+';'+str(c.ra.deg)+';'+str(c.dec.deg)+";\n")

cat=np.genfromtxt("psrcat.csv",delimiter=';',dtype=str)
cat_name=cat.T[0]
cat_p0=cat.T[1]
cat_dm=cat.T[2]
cat_ra=cat.T[3].astype(float)
cat_dec=cat.T[4].astype(float)
changedm=np.where(cat_dm=='*')
cat_dm[changedm]=0
cat_dm=cat_dm.astype(float)
########################
fname=values.fof
pos=read_fredda_hdr(values.files[0],46,47)
psrcat_format=values.column
f=open(fname,'r')
fof=f.readlines()
f.close()
foflines=np.array(fof,dtype=str)
fredda=np.loadtxt(fname)

idio=FileOperations()
idio.open_files(fname,fname)
#f1=open(fname+".psr",'w')
#f2=open(fname+".frb",'w')
for i,xy in enumerate(pos,0):
    print(i,xy)
    bm_list,fline_list=fredda_read(fredda,i)
    ra=xy[0]
    dec=xy[1]
    print("Searching "+str(len(bm_list))+" Candidates in Beam"+str(i))
    if len(bm_list) > 0:
        searcher(i,ra,dec,bm_list,foflines[fline_list],idio,cat_name,cat_dm,cat_ra,cat_dec)

#f1.close()
#f2.close()
idio.close_files()

#if __name__ == '__main__':
#    _main()
