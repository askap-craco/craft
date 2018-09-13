#!/usr/local/bin/python
import sys
import requests
import json
import logging
import socket
import datetime
import boto3

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
    parser.add_argument('--text',dest="text",type=str)
    parser.set_defaults(verbose=False)

    args = parser.parse_args()
    sb = args.sbid
    text = args.text
    sb_alias = args.alias
    scan = args.scan
    ant = args.ant
    cid = args.cid
    cid_vol = args.cid_vol
    mode = args.mode

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    
    slack (text,sb,scan,cid,ant,sb_alias)

def fix_beam(beam):
    if (beam < 10):
            beam = '0'+str(beam)
    return beam

def slack(text,sb,scan,cid,ant,sb_alias):
    
    plot_dir = '/data/TETHYS_1/bha03n/test/auto_plots'

    data = sys.stdin.readlines()
    print type(data)
    print data[0]
    print "The contents of snoopy.log are", data[1]
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
    ASKAP_SLACK_URL = "https://hooks.slack.com/services/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"
    SLACK_URL = { "askap": ASKAP_SLACK_URL }
    HOST = socket.gethostname()

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
                        }
                ]
        }
    jmessage = json.dumps(message)
    r = requests.post(ASKAP_SLACK_URL, jmessage, headers = {'content-type': 'application/json'})
    
    

if __name__ == '__main__':
        _main()
