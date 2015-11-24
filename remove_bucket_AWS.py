#!/usr/bin/env python

''' Remove an AWS S3 bucket with data in it (added with upload_data_AWS.py)
'''

from upload_download_s3 import remove_s3_bucket
from boto.s3.connection import S3Connection

bucket_name='astrocompute_testbucket0'
aws_access={'aws_access_key_id': 'AKIAJCAOONJFMHT2IIRA',
 'aws_secret_access_key': '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'}
conn=conn = S3Connection(**aws_access)

remove_s3_bucket(bucket_name,conn)
