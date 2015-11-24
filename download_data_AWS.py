#!/usr/bin/env python

'''Download data from AWS S3 bucket (added to bucket with upload_data_AWS.py)
   Notes: key_name not full file name yet so just using wildcard for now.
   		  output_dir adds that directory in place you call script for now.
'''
import sys
sys.path.append('/home/ubuntu/aws_controller')
from upload_download_s3 import download_from_s3

key_name='*'
bucket_name='BUCKET_NAME'
output_dir='/home/ubuntu/data'

download_from_s3(key_name, bucket_name, conn=None,
                     aws_access={'aws_access_key_id': 'ACCESS_KEY_HERE',
 'aws_secret_access_key': 'SECRET_ACCESS_KEY_HERE'}, output_dir=output_dir)
