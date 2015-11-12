#!/usr/bin/env python

'''Download data from AWS S3 bucket (added to bucket with upload_data_AWS.py)
   Notes: key_name not full file name yet so just using wildcard for now.
   		  output_dir adds that directory in place you call script for now.
'''

from upload_download_s3 import download_from_s3

key_name='*'
bucket_name='astrocompute_testbucket0'
output_dir='datadownload'

download_from_s3(key_name, bucket_name, conn=None,
                     aws_access={'aws_access_key_id': 'AKIAJCAOONJFMHT2IIRA',
 'aws_secret_access_key': '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'}, output_dir=output_dir)
