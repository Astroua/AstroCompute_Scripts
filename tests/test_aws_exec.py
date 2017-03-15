
'''
Test AWS execution
'''

import time
import os
import sys
import inspect
from datetime import datetime
import traceback as tr

# Append the above path for imports
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from aws_controller.upload_download_s3 import upload_to_s3, download_from_s3
from aws_controller.utils import timestring, human_time

from aws_exec import run
from utils import convert_param_format

from boto import Config

data_path = sys.argv[1]

# Make a new output folder
tstamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
output_folder = os.path.join(data_path, "output_{}".format(tstamp))

print("Starting at: " + human_time())
try:
    config = Config()
    config.load_credential_file(os.path.join(os.path.expanduser("~"),
                                             ".aws/credentials"))
    info = config.items("default")[2:]
    key = info[0][1]
    secret = info[1][1]

    aws_params = {'key': key, 'secret': secret, 'region': 'us-west-2',
                  'image_id': 'ami-b7fc75d7',
                  # use ami-cb6e75aa to pull scripts from e-koch's fork
                  'security_groups': 'launch-wizard-1',
                  'instance_type': 'm3.large'}

    params = convert_param_format("param.txt", to="dict")

    # Add in the user parameters.
    params["user_email"] = "koch.eric.w@gmail.com"
    params["user_name"] = "Dill Pickle"
    params["time_limit"] = 1
    params["summary_url"] = "https://server-name"
    params["project_email"] = "ualberta.astrocompute@gmail.com"
    start_time = timestring()

    print("Uploading at: " + human_time())
    upload_to_s3(params['target'].lower() + "_" + start_time,
                 os.path.join(data_path, 'testconcat.ms'),
                 key_prefix="data/", create_bucket=True)

    time.sleep(10)

    print("Running process at: " + human_time())
    run(start_time, params, aws_params)

    time.sleep(10)

    print("Downloading results at: " + human_time())
    download_from_s3("data_products/*",
                     params['target'].lower() + "_" + start_time,
                     output_dir=output_folder)

except Exception as e:
    print("Failed at " + human_time())
    print(e)
    tr.print_exc()
