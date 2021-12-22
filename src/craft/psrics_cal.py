#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
import os
import sys
import logging
import subprocess as sbpro
from astropy.time import Time
from pylab import *
#import fred_id
#from influxdb import InfluxDBClient

#hdr=['sn', 'tsamp', 'time', 'boxcar', 'idt', 'dm', 'beamno',
'sampno_start','sampno_end', 'idt_start', 'idt_end', 'ncands']
hdr_ics=['sn', 'tsamp', 'time', 'boxcar', 'idt', 'dm', 'beamno','mjd','sampno_start','sampno_end', 'idt_start', 'idt_end', 'ncands'] ### final header

# S/N, sampno, secs from file start, boxcar, idt, dm, beamno,mjd, sampno_start, sampno_end, idt_start, idt_end, ncands

###mock version ics hdr for testing
toff=31013-27428 ### conversion time bug??
#toff=162
#### this is for specific psrcal of vela

def loadcand(fin):
    print ('Loading', fin)
    icscand = np.loadtxt(fin, usecols=(0,1,2,3,4,5,6))
    icsmask = (icscand[:,0] >= 7) & (65 < icscand[:,5]) & (icscand[:, 5] < 70) & (icscand[:, 6] == 15)
    icscand = icscand[icsmask, :]
    idxs = np.argsort(icscand[:, 1]) # sort bty sample number
    icscand = icscand[idxs, :]
    return icscand

def pandas_read(file,hdr):
    reader=pd.read_csv(file,sep=' ',comment="#",names=hdr)
    return (reader)

def read_fredda_hdr(f):
    #x=47
    #y=48
    infile=open(f,'r')
    hdr=infile.readlines()
    infile.close()
    for i in range(len(hdr)):
        if hdr[i][0:7]=='BEAM_RA':
            x=i
        if hdr[i][0:7]=='BEAM_DE':
            y=i
        if i[0:5]=="ANTNO":
            antname="AK{0:02}".format(int(i.split(" ")[1].split("\n")[0]))
    ra=np.fromstring(hdr[x].split()[1],sep=',')[0:36]
    dec=np.fromstring(hdr[y].split()[1],sep=',')[0:36]
    pos=np.array([ra,dec]).T
    return pos,ant

def pulsar_info(psrname):
    catpsr=("psrcat -c 'dm raj decj p0' -o short_csv -nonumber -nohead "+psrname)
    catinfo=sbpro.check_output(catpsr,shell=True)
    return catinfo.decode().split(";")[:-1]

def psr_mask(antlist,psr_dm,wd=2,bm=15,dev=20):
    #######mask for pandas
    msk=(antlist.dm<psr_dm+dev)&(antlist.dm>psr_dm-dev)&(antlist.boxcar<wd)&(antlist.beamno==bm)
    return msk


def candjoin(c1, c2, t_tol=4,sntol=100):
    all_cands = []
    #ncands, ncols = c1.shape
    #print(ncands)
    for icand in c1.index:
        thecand = c1.loc[icand]  ## select line
        #print(thecand.dm)
        candt = thecand.tsamp ### tsamp value
        nearest_idx = np.argmin(abs(c2.tsamp - candt))
        nearest_cand = c2.loc[nearest_idx]
        tdist = abs(nearest_cand.tsamp - candt)
        if tdist < t_tol and thecand.sn < sntol:
            #print("yes")
            #print("\n--------------\n")
            #print (icand, candt, thecand)
            #print(nearest_idx, tdist, nearest_cand)
            #print (tdist)
            #print (thecand[1], nearest_cand[1])
            #print("found")
            all_cands.append((thecand, nearest_cand))

            #print(thecand, nearest_cand)

    return np.array(all_cands)


def candjoin_mjd(c1, c2, t_tol=(10e-3/3600),sntol=100):
    all_cands = []
    #ncands, ncols = c1.shape
    #print(ncands)
    for icand in c1.index:
        thecand = c1.loc[icand]  ## select line
        #print(thecand.dm)
        candt = thecand.mjd ### tsamp value
        nearest_idx = np.argmin(abs(c2.mjd - candt))
        nearest_cand = c2.loc[nearest_idx]
        tdist = abs(nearest_cand.mjd - candt)*3600
        if tdist < t_tol and thecand.sn < sntol:
            print("yes")
            #print("\n--------------\n")
            #print (icand, candt, thecand)
            #print(nearest_idx, tdist, nearest_cand)
            #print (tdist)
            #print (thecand[1], nearest_cand[1])
            #print("found")
            all_cands.append((thecand, nearest_cand))

            #print(thecand, nearest_cand)

    return np.array(all_cands)


#def _main():
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

if __name__ == '__main__':
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--show', action='store_true', help='Show plots')
    parser.add_argument('-g','--gen', action='store_true', help='Generate plots')
    parser.add_argument("--influx", action='store_true', help='print influx dict')
    parser.add_argument('-i','--ics',type=str,default='fredda.cand.fof',help='ics file here')
    parser.add_argument('-b','--beam',type=int,default=15,help='beam position')
    parser.add_argument("-p",'--pulsar',type=str,default='B0833-45',help='Pulsar Name for psrcat check')
    parser.add_argument("-w",'--width',type=float,default=20,help="boxcar limit")
    parser.add_argument("-m",'--mjd',type=float,default=0,help="start time")
    parser.add_argument("-x",'--sncut',type=float,default=100,help="max s/n for ics pulses")
    parser.add_argument("-t",'--tsamp',type=float,default=0.00126646875,help='time sample of filterbank (seconds)')
    parser.add_argument(dest='files', nargs='+',help='list of individual antennas')
    parser.add_argument('--mock_test1', action='store_true', help='Harry 1st dataset mocktests')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if values.mock_test1:
        hdr_ics=['sn', 'tsamp', 'time', 'boxcar', 'idt', 'dm', 'beamno','mjd','sampno_end', 'idt_start', 'idt_end', 'ncands']
        toff=31013-27428
icscand=pandas_read(values.ics,hdr_ics)
icscand.tsamp-=toff  ### time offset for this dataset
tsamp=values.tsamp
#antenna_fof=values.files
mjd=Time(values.mjd,format="mjd")
unix=mjd.unix
bm=values.beam
wd=values.width
icsmax=values.sncut
'''
antnames = [f.split('_')[2] for f in antenna_fof]
nant=len(antnames)
antcands=[pandas_read(f,hdr) for f in antenna_fof]
antnames.append('ICS')
antcands.append(icscand)
'''


antenna_fof=[f.split('/') for f in values.files]
antnames = [f[2] for f in antenna_fof]
nant=len(antnames)
antcands=[pandas_read(f,hdr_ics) for f in values.files] ##+"offline.cand.fof"
antnames.append('ICS')
antcands.append(icscand)


#pos=read_fredda_hdr(values.header)
psr_name=values.pulsar
pdm,pra,pdec,pp0=pulsar_info(psr_name)

pdm,pp0=float(pdm),float(pp0)

mask=[]
for i in antcands[:-1]:
    #print i.dm > pdm
    mask.append(psr_mask(i,pdm,wd,0,dev=20))

mask.append(psr_mask(icscand,pdm,wd,bm,dev=20))
#### stats for

    ###each antenna time-candidate
for i,ant in enumerate(antnames[:]):
    tlength=antcands[i].tsamp[mask[i]].max()
    ntol=tlength*tsamp/pp0
    print(ntol)
    print(i,ant)
    p_retr=mask[i].sum()/ntol
    can_rate=float(len(mask[i]))/(tlength*tsamp)
    if values.gen:
        plt.figure(i)
        plt.scatter(antcands[i].tsamp[mask[i]],antcands[i].sn[mask[i]])
        plt.plot(icscand.tsamp[mask[-1]],icscand.sn[mask[-1]],c='orange')
        if values.show:
            plt.show()
        else:
            plt.savefig(ant+"_temp.png")
        plt.close()
reslist=[]
#### cross match results
'''
for i in range(len(antcands[:-1])):
    ant=antnames[i]
    print("calculating ics efficiency for "+ant)
    jcands=candjoin(icscand,antcands[i],t_tol=4,sntol=icsmax)
    bb=psr_mask(jcands[:,1].T,pdm)
    icssn = jcands[:, 0, 0][bb]
    singlesn = jcands[:, 1, 0][bb]
    '''

for i in range(len(antcands[:-1])):
    #print i//5,i%5
    ant=antnames[i]
    print("calculating ics efficiency for "+ant)
    jcands=candjoin(icscand.loc[mask[-1]],antcands[i].loc[mask[i]],t_tol=4,sntol=icsmax)
    icssn = jcands[:, 0, 0]
    singlesn = jcands[:, 1, 0]
    mjd=Time(jcands[:,0,7].min(),format="mjd")
    unix=mjd.unix
    msk = icssn < 110
    fit = np.polyfit(singlesn[msk], icssn[msk], 1)
    sen_eff=fit[0]**2
    goodness=np.sqrt(np.sum(((np.polyval(fit, singlesn)-icssn)**2)/len(icssn)))
    print ("Antenna "+ant+" "+'Fit gradient={:0.1f} = {:0.1f} antennas\n Standard Deviation={:0.1f}'.format(fit[0], sen_eff,goodness))
    if values.gen:
        #print("plotting")
        figure(0)
        subplot(3,6,i+1)
        scatter(singlesn, icssn)
        plot(singlesn, np.polyval(fit, singlesn), 'r')
        xlabel("Single Antenna S/N")
        ylabel("ICS S/N")
        title(ant)
        plt.figure(1)
        plt.title("Antenna "+ant+" "+'Fit gradient={:0.1f} = {:0.1f} antennas\n Standard Deviation={:0.1f}'.format(fit[0], sen_eff,goodness))
        plt.scatter(singlesn, icssn)
        plt.plot(singlesn, np.polyval(fit, singlesn), 'r')
        plt.xlabel("Single Antenna S/N")
        plt.ylabel("ICS S/N")
        plt.savefig("gradient_"+ant+".png")
        plt.close()
#print(icscand)

    field_dict={"candidate_rate":can_rate,
                "pulse_retrival":p_retr,
                "sens":sen_eff,
                "chi":goodness,
                }
    body = {'measurement':'psrics_stat',
            'tags':{'ant':ant,'beam':bm},
            'time':unix,
            'fields':field_dict
            }
    if values.influx:
        print(body)
    reslist.append(body)
if values.show:
    plt.show()
plt.close()



#if __name__ == '__main__':
#    _main()
