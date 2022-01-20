#!/usr/bin/env python
#from pylab import *
# import matplotlib as mpl
#import matplotlib.pyplot as plt
import numpy as np
import os
import sys
#import logging
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
#from astropy.table import QTable, Table, Column

__author__ = "CRAFT Harry Qiu <hqiu0129@uni.sydney.edu.au>"


class FileOperations:
    def __init__(self,cands):
        self.cand=np.loadtxt(cands).T
        readcand=open(cands,'r')
        self.candtxt=readcand.readlines()
        readcand.close()
        self.candheader=self.candtxt[0].split('#')[1].split('\n')[0].split(',')
        self.file1 = open('sniffy.log', 'w')
        self.file1.write("#PSRname, "+self.candtxt[0])
        self.file2 = open('fifi.log', 'w')
        self.file2.write(self.candtxt[0])
    def outputfof(self,param):
        print((self.candheader[param]))
        return self.cand[param]
    def collect_multibeam(self):
        if len(self.candtxt)>2:
            self.times=np.unique(self.cand[1])
            bgcand=np.array([])
            self.newcandtxt=[]
            for i in self.times:
                idx=np.where(self.cand[1]==i)[0]
                if len(idx) >1:
                    dmcheck=np.diff(self.cand[5][idx])
                    if np.max(dmcheck) > 5:
                        continue
                print((self.cand.T[idx]))
                primarybeam=np.argmax(self.cand[0][idx])
                print(primarybeam)
                print((self.cand.T[idx][primarybeam]))
                self.newcandtxt.append(self.candtxt[idx[primarybeam]+1])
                bgcand=np.append(bgcand,self.cand.T[idx][primarybeam])
            self.reducedcand=bgcand.reshape(-1,len(self.candheader))
            return self.newcandtxt
        else:
            self.newcandtxt=self.candtxt[1:]
            self.reducedcand=self.cand

    def getparams(self,param):
        print((self.candheader[param]))
        return self.reducedcand.T[param]
    def write_psr(self,line,psrname,prnt=False):
        self.file1.write(psrname+"\t"+self.newcandtxt[line])
        if prnt==True:
            print(line)
    def write_frb(self,line,prnt=False):
        self.file2.write(self.newcandtxt[line])
        if prnt==True:
            print(line)
    def close_files(self):
        self.file1.close()
        self.file2.close()

class PsrPnt:
    def __init__(self,coords,dm,radius=1,coord_units=u.degree,rad_unit=u.degree):
        self.coords=coords
        self.beamsize=radius*rad_unit
        self.coord_units=coord_units
        self.rad_unit=rad_unit
        self.dm=dm#*(u.parsec*(u.cm)**-3)
        MIN_FLOAT = sys.float_info[3]
        psrcat=os.environ['FRED_CAT']
        if os.path.exists(psrcat):
            cat=np.genfromtxt(psrcat,delimiter=';',dtype=str)
            self.cat_name=cat.T[0]
            self.cat_p0=cat.T[1]
            cat_dm=cat.T[2]
            cat_ra=cat.T[3].astype(float)
            cat_dec=cat.T[4].astype(float)
            self.psrcat=SkyCoord(ra=cat_ra*coord_units, dec=cat_dec*coord_units)
            changedm=np.where(cat_dm=='*')
            cat_dm[changedm]=-10
            self.cat_dm=cat_dm.astype(float)

        else:
            raise NameError ("No psrcat.csv!")

        #print(psrcat)
    def coordmatch(self,ra,dec,beamradius=1.5):
        # self.ra=ra
        # self.dec=dec
        location=SkyCoord(ra=ra*self.coord_units, dec=dec*self.coord_units)
        sep=location.separation(self.psrcat)
        mask=sep<(self.beamsize*beamradius)
        #print(self.cat_name[mask],self.psrcat[mask])
        return mask
    def dmmatch(self,dm,dmerr=1):
        #dmerr=np.abs(dmerr)
        mask= (dm-dmerr<self.cat_dm)*(self.cat_dm<dm+dmerr)
        #print(self.cat_name[mask],self.cat_dm[mask],self.psrcat[mask])
        return mask
    def full_crossmatch(self,ra,dec,dm,dmerr=1,beamradius=1.5):
        print((dm,ra,dec))
        coordmask=self.coordmatch(ra,dec,beamradius)
        dmmask=self.dmmatch(dm,dmerr)
        mask=coordmask*dmmask
        location=SkyCoord(ra=ra*self.coord_units, dec=dec*self.coord_units)
        sep=location.separation(self.psrcat[mask])
        print("PSRCAT results as follows:")
        print((self.cat_name[mask],self.cat_dm[mask],self.psrcat[mask]))
        return mask,sep
    def match_all(self,dmerr=5,beamradius=1.5):
        k=0
        printids=[]
        if isinstance(self.dm,np.float64):
            print(("Candidate "+str(k)))
            i=self.coords
            j=self.dm
            pmask,psep=self.full_crossmatch(i[0],i[1],j,dmerr,beamradius)
            if len(psep)>0:
                print("preview pulsar")
                print((self.cat_name[pmask][np.argmin(psep)]))
                printids.append([k,self.cat_name[pmask][np.argmin(psep)]])
            else:
                printids.append([k,'None'])
        else:
            for i,j in zip(self.coords,self.dm):
                print(("Candidate "+str(k)))
                pmask,psep=self.full_crossmatch(i[0],i[1],j,dmerr,beamradius)
                if len(psep)>0:
                    print("preview pulsar")
                    print((pmask[np.argmin(psep)]))
                    printids.append([k,self.cat_name[pmask][np.argmin(psep)]])
                else:
                    printids.append([k,'None'])
                k+=1
        return printids



    def readmask(self,mask):
        print((self.cat_name[mask],self.cat_dm[mask],self.psrcat[mask]))





class hdrfiles():
    def read_beampos(self):
        hdr=self.hdr
        for i in range(len(hdr)):
            if hdr[i][0:7]=='BEAM_RA':
                x=i
            if hdr[i][0:7]=='BEAM_DE':
                y=i
        #print x,y
        #sbid=hdr[9].split()[1]
        #obsid=hdr[8].split()[1]
        #ant=hdr[22].split()[1]
        ra=np.fromstring(hdr[x].split()[1],sep=',')[0:36]
        dec=np.fromstring(hdr[y].split()[1],sep=',')[0:36]
        pos=np.array([ra,dec]).T
        #self.beampos=pos
        ######/data/TETHYS_1/craftop/auto_cands//SB01231_20180329103213_co18_tethys4.cand.fof.mbeam.2.good
        ####this isn't finished need fof file name format generation "SB{0:05}_".format(int(sbid))+" _co18_"+".cand.fof.mbeam.2.good"
        return pos
    def __init__(self,filename):
        infile=open(filename,'r')
        self.hdr=infile.readlines()
        self.beampos=self.read_beampos()
        infile.close()
    def get_beampos(self,beamno):
        return(self.beampos[beamno])

### TEST run below
##hdrf,freddafof='ICS_SB13806_C000.hdr','test2.cand.fof'
##dmerr,beamradius,beam=10,5,1
def freddachecker(hdrf,freddafof,dmerr=10,beamradius=5,beam=1):
    candidates=FileOperations(freddafof)
    print((freddafof,hdrf))
    candidates.collect_multibeam()
    cand_beams=candidates.getparams(6).astype(np.int)
    cand_dms=candidates.getparams(5)
    header=hdrfiles(hdrf)
    coords=header.get_beampos(cand_beams)
    print(coords)
    dmlist=cand_dms
    print(dmlist)
    sniffer=PsrPnt(coords,dmlist,radius=beam)
    outlist=sniffer.match_all(dmerr,beamradius)
    print("outlist produced")
    print(outlist)
    for i in outlist:
        # print(idx)
        idx=i[0]
        print(idx)
        psrs=i[1]
        print(psrs)
        if psrs=='None':
            candidates.write_frb(idx,prnt=True)
        else:
            candidates.write_psr(idx,psrs,prnt=True)
    candidates.close_files()

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-r','--radius',type=float,default=1.5,help='search radius')
    #parser.add_argument('-x','--sncut',type=float,default=10,help='fredda.frb file threshold for sn')
    parser.add_argument('-d','--dmlim',type=float,default=5,help='dm errorbar')
    parser.add_argument('-f','--candlist',type=str,default='fredda.cand.fof',help='fof file here') #### fredda candidate file input here
    parser.add_argument(dest='files', nargs='+',help='hdrfile')
    values = parser.parse_args()
    hdrf=values.files[0]
    freddafof=values.candlist
    freddachecker(hdrf,freddafof,dmerr=values.dmlim,beamradius=values.radius,beam=0.9)

if __name__ == '__main__':
    _main()
