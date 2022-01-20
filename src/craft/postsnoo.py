#!/usr/bin/env python

import sys
import requests
import json
import logging
import socket
import datetime
import boto3
import os

ACCESS_KEY = os.getenv('CRAFT_S3_ACCESS')
SECRET_KEY = os.getenv('CRAFT_S3_SECRET')
BUCKET = os.getenv('CRAFT_S3_BUCKET')
ASKAP_SLACK_URL = os.getenv('ASKAP_SLACK_URL')

def _main ():
    
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('--sbid', dest="sbid", type=str)
    parser.add_argument('--scan', dest="scan", type=str)
    parser.add_argument('--alias',dest="alias",type=str)
    parser.add_argument('--ant',dest="ant",type=str)
    parser.add_argument('--cid',dest="cid",type=str)
    parser.add_argument('--cid_vol',dest="cid_vol",type=str)
    parser.add_argument('--mode',dest="mode",type=int)
    parser.add_argument("--attachment", type=str, help="Full path of the image")
    parser.add_argument('--text',dest="text",type=str)
    parser.set_defaults(verbose=False)

    values = parser.parse_args()
    sb = values.sbid
    text = values.text
    alias = values.alias
    scan = values.scan
    ant = values.ant
    cid = values.cid
    attachment = values.attachment
    #cid_vol = args.cid_vol
    #mode = args.mode

    # Temp fix
    #sb = sb[0:2]+'0'+sb[2:]
    cid = fix_cid(cid)
        
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    #slack(values,sb,cid,text)
    slack (text,sb,scan,cid,ant,alias,attachment)

def fix_cid(cid):
    if (len(cid) == 1):
        cid = 'C00'+cid
    elif (len(cid) ==2):
        cid = 'C0'+cid
    else:
        cid = 'C'+cid
    return cid

def fix_beam(beam):
    if (beam < 10):
            beam = '0'+str(beam)
    return beam

def slack(text,sb,scan,cid,ant,alias,attachment):
    
    
    plot_dir = '/data/TETHYS_1/bha03n/test/auto_plots'

    data = sys.stdin.readlines()
    print(type(data))
    print(data[0])
    print("The contents of snoopy.log are", data[1])
    line = data[1]
    values = line.split(" ")
    snr = values[0]
    sampno = values[1]
    t_start = values[2]
    width = values[3]
    dm = values[5]
    beam = int(values[6])
    mjd = values[7]
    latency_ms = values[8]

    beam = fix_beam(beam)
    
    #image1 = str(sb) + "_" + str(scan) + "_" + str(cid) + "_" + str(ant) + '_' + str(beam) + '_on.png'
    #image2 = str(sb) + "_" + str(scan) + "_" + str(cid) + "_" + str(ant) + '_' + str(beam) + '_off.png'
    
    SLACK_URL = { "askap": ASKAP_SLACK_URL }
    HOST = socket.gethostname()
   
    
    now = datetime.datetime.now()
    expires = now + datetime.timedelta(minutes=1000080)

    if (attachment != None):
        image_key = attachment
        print("Supplied attachment image is",image_key)
        client = boto3.client('s3',aws_access_key_id=ACCESS_KEY,aws_secret_access_key=SECRET_KEY)
        client.upload_file(image_key, BUCKET, image_key,ExtraArgs={'ContentType': "image/png", 'ACL': "public-read",'Expires': expires})
        url = client.generate_presigned_url(ClientMethod='get_object',ExpiresIn=604800,Params={'Bucket':BUCKET,'Key':image_key})

    else:
        print("No attachement is provided")
        url = None
    '''    
    client = boto3.client('s3',aws_access_key_id=ACCESS_KEY,aws_secret_access_key=SECRET_KEY)
    
    image_key1 = plot_dir + "/" + image1
    client.upload_file(image_key1, BUCKET, image_key1,ExtraArgs={'ContentType': "image/png", 'ACL': "public-read",'Expires': expires})
    url1 = client.generate_presigned_url(ClientMethod='get_object',Params={'Bucket':BUCKET,'Key':image_key1})
    '''
    

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
            				"value": alias,			
                                        "short": True
                                    },

                                    {
                                        "title":"Scan ID",
                                        "value": scan,
                                        "short": True
                                    },
                                    {
                                        "title": "Capture ID",
                                        "value": str(cid),
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
                                        "title":"Latency (ms)",
                                        "value": latency_ms,
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
