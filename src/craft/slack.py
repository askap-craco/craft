#!/usr/bin/env python
"""
Scripts for posting stuff on slack

Copyright (C) CSIRO 2015
"""

import sys
import requests
import json
import logging
import socket
import datetime
import boto3
import argparse


def _main():
    from argparse import FileType, ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='fields', nargs='+')
    parser.add_argument('--text',type=str)
    parser.add_argument("--attachment", type=str, help="Full path of the image")
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    fields = []
    
    for s in values.fields:
        print(s)
        bits = s.split('=')
        k,v = bits
        f = {'title':k,
             'value':v,
             'short':True
        }

        fields.append(f)

    ACCESS_KEY = 'AKIAIUYCJ7ATHZQ274XA'
    SECRET_KEY = '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'
    BUCKET = 'casstestaskapmro2'

    ASKAP_SLACK_URL = "https://hooks.slack.com/services/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"
    SLACK_URL = { "askap": ASKAP_SLACK_URL }
    HOST = socket.gethostname()

    now = datetime.datetime.now()
    expires = now + datetime.timedelta(minutes=10080)
    
    if (values.attachment != None):
            image_key = values.attachment
            print("Supplied attachment image is",image_key)
            
            client = boto3.client('s3',aws_access_key_id=ACCESS_KEY,aws_secret_access_key=SECRET_KEY)
            client.upload_file(image_key, BUCKET, image_key,ExtraArgs={'ContentType': "image/png", 'ACL': "public-read",\
                                                                       'Expires': expires})
            url = client.generate_presigned_url(ClientMethod='get_object',ExpiresIn=604800,Params={'Bucket':BUCKET,'Key':image_key})
    else:
            print("No attachement is provided")
            url = None
            

    message = {"text":values.text, "attachments" :
               [
                   {
                       "fields": fields,
                       "image_url": url
                    }
               ]
               }
    
    jmessage = json.dumps(message)
    r = requests.post(ASKAP_SLACK_URL, jmessage, headers = {'content-type': 'application/json'})
if __name__ == '__main__':
    _main()
