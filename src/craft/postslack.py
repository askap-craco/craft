#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import os
import sys
import logging
import boto3
import json
import requests
import datetime

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

ACCESS_KEY = os.getenv('CRAFT_S3_ACCESS')
SECRET_KEY = os.getenv('CRAFT_S3_SECRET')
BUCKET = os.getenv('CRAFT_S3_BUCKET')
ASKAP_SLACK_URL = os.getenv('ASKAP_SLACK_URL')

def upload_file(fname):
    image_key = os.path.abspath(fname)
    now = datetime.datetime.now()
    expires = now + datetime.timedelta(minutes=1000080)

    client = boto3.client('s3',aws_access_key_id=ACCESS_KEY,aws_secret_access_key=SECRET_KEY)
    client.upload_file(image_key, BUCKET, image_key,ExtraArgs={'ContentType': "image/png", 'ACL': "public-read",'Expires': expires})
    url = client.generate_presigned_url(ClientMethod='get_object',ExpiresIn=604800,Params={'Bucket':BUCKET,'Key':image_key})
    return url

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-i','--image', help='Attach image')
    parser.add_argument('-t','--text', help='Add text', default='')
    parser.add_argument(dest='fields', nargs='*')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    # get text
    if not os.isatty(sys.stdin.fileno()):
        values.text += sys.stdin.read()

    fields = []
    for s in values.fields:
        bits = s.split('=')
        if len(bits) != 2:
            raise ValueError('Cannot parse field {}'.format(s))
        k,v = bits
        f = {'title':k,
             'value':v,
             'short':True
        }

        fields.append(f)
        
    image_url = None
    if values.image:
        image_url = upload_file(values.image)

    message = {}
    if len(values.text) > 0:
        message['text'] = values.text

    if image_url or len(fields) > 0:
        attachment = {}
        if len(fields) > 0:
            attachment['fields'] = fields
            
        if image_url:
            attachment['image_url'] = image_url

            
        message['attachments'] = [attachment]

    
    jmessage = json.dumps(message)
    r = requests.post(ASKAP_SLACK_URL, jmessage, headers = {'content-type': 'application/json'})

if __name__ == '__main__':
    _main()
