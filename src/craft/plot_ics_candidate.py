#!/usr/bin/env python

import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import glob
from subprocess import Popen,PIPE,STDOUT,call
import datetime

# for slack
import requests
import json
import socket
import boto3

def _main():

    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--sbid')
    parser.add_argument('--scan')
    parser.add_argument('--alias')
    parser.add_argument('--ant')
    parser.add_argument('--cid')
    parser.add_argument('--text')
    parser.add_argument('--mode', help='CRAFT mode')
    parser.add_argument('--fmode', help ='filterbank to process. Choose on/off')
    parser.set_defaults(verbose=False)

    args = parser.parse_args()
    sb = args.sbid
    sb_alias = args.alias
    scan = args.scan
    ant = args.ant
    cid = args.cid
    text = args.text
    fmode = args.fmode
                    
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    check_realcand_ics(args)
    slack (args)

    
def fix_sb(sb):
    sb_new = sb[:2] + sb [3:7]
    return sb_new

def pol_beam(beam):
    
    beam = int(beam)
    beamx = 2 * beam
    beamy = 2 * beam + 1
    
    if (beamx < 9):
        beamx = '0'+str(beamx)
                    
    if (beamy < 9):
        beamy = '0'+str(beamy)

    return str(beamx),str(beamy)


def fix_beam(beam):

    if (beam < 10):
        beam = '0'+str(beam)
    return str(beam)
    
def make_vcraft_plots(sb,scan,cid,cid_vol,dm,beam_v,beam,ant,mjd,plot_time):
    # Here cid would be the capture id for the particular observation.
    # Here beam should be vcraft beam for the particular observation.
    
    print("X polarisation is",beam_v[0])
    print("Y polarisation is",beam_v[1])
    print("Making plots using filterbank made from voltages")
    print("Antenna to work on is",ant)
    print(type(ant))
    original_sb = sb
    sb = fix_sb(sb)
    plot_dir = '/data/TETHYS_1/bha03n/test/auto_plots/'
    
    for i in range (len(beam_v)):
        print("I is",i)
        fil_path = ant +'/' + str(beam_v[i]) + '/beam.fil' 
        vcraft_path = ant+ '/'  + str(beam_v[i]) +'/'
        
        my_file = os.path.isfile(fil_path)

        if my_file:
            print("file is present")
        else:
            print("Making a beam.fil")
            cmd = ' numactl --cpunodebind 0 /home/craftop/craftdev/python/vcraft2fil.py -i 4 ' + vcraft_path + '*.vcraft' + ' -o ' + vcraft_path+ 'beam.fil'
            print(cmd)
            os.popen(cmd)
         
         
        assert os.path.isfile(fil_path), 'Looks like vcraft2fil didnt work'
        
        cmd1 = 'dspsr -cepoch=start   -D ' + str(dm) + ' -c ' + str(plot_time) + '  -T ' + str(plot_time) + ' -k askap -N source -O '+ plot_dir +original_sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam_v[i]) + '_on'+' ' + fil_path
        print(cmd1)
        os.popen(cmd1)

    print("Adding both the polarisations")         
    cmd2 = 'psradd -o ' + plot_dir +original_sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam) + '_on_add.ar ' +  plot_dir +original_sb+'_'+scan+'_'+str(cid)+'_'+ant+'_beam*_on.ar'
    print(cmd2)
    os.popen(cmd2)
    
    print("Now plotting")
    cmd3 = 'psrplot -p freq+ -c psd=0 -c x:unit=ms -jD -j "F 54"  -D ' + plot_dir + original_sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam) + '_on'+'.png/png ' + plot_dir + original_sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam) + '_on_add.ar'
    print(cmd3)
    os.popen(cmd3)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def incoherent_sum(sb,scan,cid,beam,flag):
    
    plot_dir = '/data/TETHYS_1/bha03n/test/auto_plots/'
    if (flag==1):
    
        print("Using filterbanks created from voltages") 
        

        in_file_on = '{}_{}_{}_*{}_on_add.ar'.format(sb, scan, cid, beam)
        out_file_on = '{}_{}_{}_{}_on_ics.ar'.format(sb, scan, cid, beam)
        plot1_on = '{}_{}_{}_{}_on_ics.DFTp'.format(sb, scan, cid, beam)

        inf_on = os.path.join(plot_dir,in_file_on)
        outf_on = os.path.join(plot_dir,out_file_on)
        plot1f_on = os.path.join(plot_dir,plot1_on)
        cmd = 'psradd  -o ' + outf_on + ' ' + inf_on
        print(cmd)
        os.popen(cmd)
        cmd2 = 'pav -DFTp -N 3,3  ' + inf_on + ' ' + outf_on + ' -g ' + plot1f_on + '.png/PNG'
        print(cmd2)
        os.popen(cmd2)
    else:
        
        print(" Using offline filterbanks")
        in_file_off = '{}_{}_{}_*{}_off.ar'.format(sb, scan, cid, beam)
        out_file_off = '{}_{}_{}_{}_off_ics.ar'.format(sb, scan, cid, beam)
        plot1_off = '{}_{}_{}_{}_off_ics.DFTp'.format(sb, scan, cid, beam)

        inf_off = os.path.join(plot_dir,in_file_off)
        outf_off = os.path.join(plot_dir,out_file_off)
        plot1f_off = os.path.join(plot_dir,plot1_off)
        cmd3 = 'psradd  -o ' + outf_off + ' ' + inf_off
        print(cmd3)
        os.popen(cmd3)
        cmd4 = 'pav -DFTp -N 3,3  ' + inf_off + ' ' + outf_off + ' -g ' + plot1f_off + '.png/PNG'
        print(cmd4)
        os.popen(cmd4)
    
def make_plots(sb,scan,cid,dm,beam,ant,mjd,plot_time):
	# Here CID should be the scan ID and beam should be the actual beam

    print("Making plots using offline filterbank")
    fil_path = '/data/TETHYS_1/craftop/data/'+ sb + '/' + scan + '/' + ant + '/' + str(cid) + '/*' + str(beam) + '.fil'
    
    plot_dir = '/data/TETHYS_1/bha03n/test/auto_plots/'
    start_pulse,nsamp,tsamp = find_start_time(sb,scan,cid,beam,mjd,ant)
    print(nsamp)
    print(tsamp)
    tsamp_s = tsamp*1e-6
    print(tsamp_s)
    tobs = nsamp*tsamp_s
    t_left = tobs - start_pulse
    
    band = 336
    freq = 1.4
    delay = 8.3 * dm * band * np.power(freq,-3)
    delay_s = delay*1e-6
    
    
    if (plot_time > t_left):
         plot_time = t_left
         
    timeoff = start_pulse - (plot_time/2)
    
    
    cmd1 = 'dspsr -cepoch=start -S ' + str(timeoff) + ' -D ' + str(dm) + ' -c ' + str(plot_time) + '  -T ' + str(plot_time) + ' -k askap -N source -O '+ plot_dir +sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam)+'_off' + ' ' + fil_path
    print(cmd1)
    cmd2 = 'psrplot -p freq+ -c psd=0 -c x:unit=ms -jD   -D ' + plot_dir + sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam)+'_off' + '.png/png ' + plot_dir + sb+'_'+scan+'_'+str(cid)+'_'+ant+'_'+str(beam) + '_off.ar'
    print(cmd2)
    os.popen(cmd1)
    os.popen(cmd2)

def check_realcand_ics(args):
    
    sb = args.sbid
    sb_alias = args.alias
    scan = args.scan
    cid = args.cid
    text = args.text
    fmode = args.fmode
    
    # I should be inside /data/TETHYS_1/craftop/data/testing/voltages/SB1783/20180605075255/ICS/C0
    ant = []
    beam_v = []
    print("Snoopy real-time candidate")
    cmd = 'pwd | sed  "s/\// /g" | awk "{print $10}"'
    capture = os.popen(cmd).read()
    capture_id = capture[-3:]
    print("CID is", cid)
    
    beam_path = 'co*/beam*'
    print(beam_path)
    for b in glob.glob(beam_path):
        print("File is",b)
        ant.append(b.split('/')[0])
        beam_v.append(b.split('/')[1])
        
    ant = np.unique(ant)
    beam_v = np.unique(beam_v)
    print("Antennas are",ant)
    print("The two polarisation beams are",beam_v)
    snop = np.loadtxt("snoopy.log")
    snr = snop[0]
    fredsnop_mjd = snop[7]
    beam = int(snop[6])
    dm = snop[5]

    mode_plot_time_map= {"0":0.7, "1":1.4,"2":3,"3":14}
    plot_time = mode_plot_time_map[args.mode]

    if (fmode=='on'):
        for i in range(len(ant)):
            print("Making plots using online filterbanks now for ",ant[i])
            flag = 1
            make_vcraft_plots(sb,scan,cid,capture_id,dm,beam_v,beam,str(ant[i]),fredsnop_mjd, plot_time)
        incoherent_sum(sb,scan,cid,beam,flag)

    else:
        for i in range(len(ant)):
            print("Making plots using offline filterbanks now for ",ant[i])
            flag = 0
            make_plots(sb,scan,cid,dm,beam,str(ant[i]),fredsnop_mjd, plot_time)
        incoherent_sum(sb,scan,cid,beam,flag)
    
def find_start_time(sb,scan,cid,beam,mjd,ant):
    # Get the MJDs from the filterbank files.
    fil_path = '/data/TETHYS_1/craftop/data/'+ sb + '/' + scan + '/' + ant + '/' + cid + '/*' + str(beam) + '.fil'  
    print("filtebank file is",fil_path)
    header_path = '/home/sha355/bin/header'
    cmd1 = header_path + ' ' + fil_path + ' -tstart'
    fil_mjd = os.popen(cmd1).read()
    print("Fil MJD is",fil_mjd)
    fil_mjd = float(fil_mjd)
    cmd2 = header_path + ' ' + fil_path + ' -nsamples'
    nsamp = os.popen(cmd2).read()
    nsamp = int(nsamp)

    cmd3 = header_path + ' ' + fil_path + ' -tsamp'
    tsamp = os.popen(cmd3).read()
    tsamp = float(tsamp)

    start_time = (mjd-fil_mjd)*86400.
    print("Start time of the pulse is",start_time)
    return start_time,nsamp,tsamp

def slack(args):
    sb = args.sbid
    sb_alias = args.alias
    scan = args.scan
    ant = 'ALL'
    cid = args.cid
    text = args.text
    fmode = args.fmode
    plot_dir = '/data/TETHYS_1/bha03n/test/auto_plots'

    snop = np.loadtxt("snoopy.log")
    snr = snop[0]
    sampno = snop[1]
    t_start = snop[2]
    width = snop[3]
    dm = snop[5]
    beam = int(snop[6])
    mjd = snop[7]
    latency_ms = snop[8]
    
    beam = fix_beam(beam)
    print("Fmode is",fmode)
    if (fmode == 'on'):
        image = '{}_{}_{}_{}_on_ics.DFTp.png'.format(sb, scan, cid, beam)
    else:
        image = '{}_{}_{}_{}_off_ics.DFTp.png'.format(sb, scan, cid, beam)

    ACCESS_KEY = 'AKIAIUYCJ7ATHZQ274XA'
    SECRET_KEY = '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'
    BUCKET = 'casstestaskapmro2'
    
    ASKAP_SLACK_URL = "https://hooks.slack.com/services/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"
    SLACK_URL = { "askap": ASKAP_SLACK_URL }
    HOST = socket.gethostname()
    
    
    now = datetime.datetime.now()
    expires = now + datetime.timedelta(minutes=10080)
    
    client = boto3.client('s3',aws_access_key_id=ACCESS_KEY,aws_secret_access_key=SECRET_KEY)
    
    image_key = plot_dir + "/" + image
    print(image_key)
    
    client.upload_file(image_key, BUCKET, image_key,ExtraArgs={'ContentType': "image/png", 'ACL': "public-read",'Expires': expires})
    url = client.generate_presigned_url(ClientMethod='get_object',Params={'Bucket':BUCKET,'Key':image_key})
    
    message = {"text": text , "attachments" :
        [
         {
         "fields": [
                    {
                    "title": "SBID",
                    "value": str(sb),
                    "short": True
                    },
                    
                    {
                    "title" :"alias",
                    "value": str(sb_alias),
                    "short": True
                    },
                    
                    {
                    "title":"Scan ID",
                    "value": str(scan),
                    "short": True
                    },
                    
                    {
                    "title":"Antenna",
                    "value": str(ant),
                    "short": True
                    },
                    
                    {
                    "title":"Beam",
                    "value": str(beam),
                    "short": True
                    },
                    
                    {
                    "title":"S/N",
                    "value": snr,
                    "short": True
                    },
                    
                    {
                    "title":"DM",
                    "value": dm,
                    "short": True
                    },
                    {
                    "title":"width",
                    "value": width,
                    "short": True
                    },
                    
                    
                    ],
         "image_url": url
         }
         ]
    }
    jmessage = json.dumps(message)
    r = requests.post(ASKAP_SLACK_URL, jmessage, headers = {'content-type': 'application/json'})
    '''
    image_key2 = plot_dir + "/" + image2
    client.upload_file(image_key2, BUCKET, image_key2,ExtraArgs={'ContentType': "image/png", 'ACL': "public-read",'Expires': expires})
    url2 = client.generate_presigned_url(ClientMethod='get_object',Params={'Bucket':BUCKET,'Key':image_key2})
    
    message2 = {"text":"Plot made using the offline filterbank", "attachments": [ {"image_url": url2 }] }
    jmessage2 = json.dumps(message2)
    r = requests.post(ASKAP_SLACK_URL, jmessage2,headers={'content-type': 'application/json'})
    '''

if __name__ == '__main__':
        _main()
    
