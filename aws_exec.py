
import json
from time import sleep

import boto

from aws_controller.controller import WORKER_SCRIPT
from aws_controller.utils import timestring
from aws_controller.launch_instance import launch


def run(timestamp, param_file, aws_file, s3_bucket=None, log_file=None):
    '''
    Run Timing scripts on AWS. See README for more info.

    A parameters file is passed into the function the controls the stages of
    the processing. The script then uploads the data

    This function is intended to be used by itself, or through a web server
    where execution requests are submitted.

    Parameters
    ----------
    timestamp : str
        Timestamp created when the request was submitted.
    param_file : str
        JSON file with parameter settings.
    aws_file : str
        JSON file with settings for the AWS execution.
    '''

    # Read the param_file file
    with open(param_file) as f:
        params = json.load(f)

    # Read the AWS settings file.
    with open(aws_file) as f:
        aws_settings = json.load(f)

    proc_name = params['target'] + "_" + timestamp
    key = aws_settings['key']
    secret = aws_settings['secret']
    region = aws_settings['region']

    # Create the queue and submit the message.
    queue_name = proc_name
    resp_name = "resp_" + proc_name

    queue = boto.sqs.connect_to_region(region,
                                       aws_access_key_id=key,
                                       aws_secret_access_key=secret).create_queue(queue_name)

    msg = json_message(params, aws_settings, proc_name)
    queue.write(queue.new_message(body=msg))

    # Launch instance

    user_data = WORKER_SCRIPT % \
        {"USER": 'ubuntu',
         "QUEUE_NAME": proc_name,
         "KEY": key,
         "SECRET": secret,
         "REGION": region,
         "RESP_QUEUE_NAME": resp_name}

    inst = launch(key, region=region, image_id=aws_settings['image_id'],
                  instance_type=aws_settings['instance_type'],
                  security_groups=aws_settings['security_groups'],
                  initial_check=False, user_data=user_data)

    # Wait for return message
    while inst.state == "running":
        sleep(60)

    # Return success/failure


def json_message(params, aws_settings, proc_name, s3_bucket=None):
    '''
    Create the json message to feed to the worker instance.
    '''

    # Command to run
    cmd = ["/usr/local/bin/CASA/casa-release-4.3.1-el6/casa", "--nologger",
           "--logfile", proc_name+".log",
           "/home/ubuntu/AstroCompute_Scripts/casa_timing_script.py",
           "/home/ubuntu/data/params.txt", "/home/ubuntu/"]

    mess = {}

    mess['proc_name'] = proc_name
    mess['key_name'] = "data/*"

    if s3_bucket is not None:
        mess['bucket'] = s3_bucket
    else:
        mess['bucket'] = proc_name

    mess["command"] = cmd

    return json.dumps(mess)


if __name__ == "__main__":
    pass
