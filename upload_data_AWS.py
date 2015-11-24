#!/usr/bin/env python

''' Upload data to AWS S3 bucket (which can then be downloaded with download_data_AWS.py)
	Notes: No capitals allowed in bucket name
'''
import sys
sys.path.append('PATH_TO_AWS_CONTROLLER')
from upload_download_s3 import upload_to_s3

bucket_name='astrocompute_testbucket0'
upload_item='/projects/mnt/JACPOT/JACPOT_V404_Cyg/VLA/lustre/aoc/ftp/e2earchive/stage/QK3610/v404_jun22_A_Cc5.ms'

upload_to_s3(bucket_name, upload_item,
                 create_bucket=True, chunk_size=52428800, conn=None,
                 aws_access={'aws_access_key_id': 'ACCESS_KEY_HERE',
 'aws_secret_access_key': 'SECRET_ACCESS_KEY_HERE'}, replace=False)
