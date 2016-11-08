
'''
Test AWS execution
'''

import time
import datetime
import os
import traceback as tr

from aws_controller.upload_download_s3 import upload_to_s3, download_from_s3
from aws_controller.utils import timestring, human_time

from aws_exec import run
from utils import convert_param_format

from boto import Config


print("Starting at: " + human_time())
try:
    config = Config()
    config.load_credential_file(os.path.join(os.path.expanduser("~"),
                                             ".aws/credentials"))
    info = config.items("default")[2:]
    key = info[0][1]
    secret = info[1][1]

    aws_params = {'key': key, 'secret': secret, 'region': 'us-west-2',
                  'image_id': 'ami-bd0219dc',
                  # use ami-cb6e75aa to pull scripts from e-koch's fork
                  'security_groups': 'launch-wizard-1',
                  'instance_type': 'm3.large'}

    params = convert_param_format("param.txt", to="dict")

    start_time = timestring()

    print("Uploading at: " + human_time())
    upload_to_s3(params['target'].lower() + "_" + start_time,
                 '/Users/atetarenk/Desktop/AstroCompute_Nov/v404_jun22_B_Cc7_bp.ms',
                 key_prefix="data/", create_bucket=True)

    time.sleep(10)

    print("Running process at: " + human_time())
    run(start_time, params, aws_params)

    time.sleep(10)

    print("Downloading results at: " + human_time())
    download_from_s3("data_products/*",
                     params['target'].lower() + "_" + start_time,
                     output_dir="/Users/atetarenk/Desktop/AstroCompute_Nov/")

except Exception as e:
    print("Failed at " + human_time())
    print(e)
    tr.print_exc()
