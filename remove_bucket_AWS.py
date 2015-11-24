#!/usr/bin/env python

''' Remove an AWS S3 bucket with data in it (added with upload_data_AWS.py)
'''
import sys
sys.path.append('PATH_TO_AWS_CONTROLLER')
from upload_download_s3 import remove_s3_bucket
from boto.s3.connection import S3Connection

bucket_name='BUCKET_NAME'
aws_access={'aws_access_key_id': 'ACCESS_KEY_HERE',
 'aws_secret_access_key': 'SECRET_ACCESS_KEY_HERE'}
conn=conn = S3Connection(**aws_access)

remove_s3_bucket(bucket_name,conn)
