
import json
from time import sleep, time

import boto
from boto import sqs

from aws_controller.controller import WORKER_SCRIPT
from aws_controller.utils import timestring
from aws_controller.launch_instance import launch


def run(timestamp, param_file, aws_file):
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
    param_file : str or dict
        JSON file or a dictionary with parameter settings.
    aws_file : str or dict
        JSON file or a dictionary with settings for the AWS execution.
    '''

    # Read the param_file file
    if isinstance(param_file, dict):
        params = param_file
    else:
        with open(param_file) as f:
            params = json.load(f)

    # Read the AWS settings file.
    if isinstance(aws_file, dict):
        aws_settings = aws_file
    else:
        with open(aws_file) as f:
            aws_settings = json.load(f)

    proc_name = params['target'] + "_" + timestamp
    key = aws_settings['key']
    secret = aws_settings['secret']
    region = aws_settings['region']

    # Create the queue and submit the message.
    queue_name = proc_name
    resp_queue_name = "resp_" + proc_name

    queue = sqs.connect_to_region(region,
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
         "RESP_QUEUE_NAME": resp_queue_name}

    inst = launch(key, region=region, image_id=aws_settings['image_id'],
                  instance_type=aws_settings['instance_type'],
                  security_groups=aws_settings['security_groups'],
                  initial_check=False, user_data=user_data)

    # Wait for return message
    # Total wait time set to 1 day
    time_max = 3600 * 24
    time_wait = 600
    t0 = time()
    while time() < t0 + time_max:
        update = inst.update()
        if update in [u"stopping", u"stopped"]:
            break
        sleep(time_wait)
    else:
        print("Reached time limit. Terminating.")

    inst.terminate()

    # Connect to the response queue
    resp_queue = sqs.connect_to_region(region,
                                       aws_access_key_id=key,
                                       aws_secret_access_key=secret).create_queue(resp_queue_name)

    # Check the response message
    mess = resp_queue.read(10)
    if mess is not None:
        content = json.loads(mess.get_body())
        print("Saving content.")
        with open("tests/test_response.txt", "w") as f:
            json.dump(content, f)

    # Clean-up queues
    queue.delete()
    resp_queue.delete()

    # Return success/failure


def json_message(params, aws_settings, proc_name):
    '''
    Create the json message to feed to the worker instance.
    '''

    # Commands to run
    clone_cmd = \
        'su - ubuntu -c "cd /home/ubuntu/; " \
        "/usr/bin/git clone git@github.com:Astroua/AstroCompute_Scripts.git"'
    timing_cmd = \
        "/usr/local/bin/CASA/casa-release-4.3.1-el6/casa --nologger " \
        "--logfile " + proc_name + "_timing.log -c "\
        "/home/ubuntu/AstroCompute_Scripts/casa_timing_script.py "\
        "/home/ubuntu/data/params.txt /home/ubuntu/"

    # Check whether Object detection should be run
    if params['runObj'] == "T":
        clean_cmd = \
            "/usr/local/bin/CASA/casa-release-4.3.1-el6/casa --nologger " \
            "--logfile " + proc_name + "_timing.log -c "\
            "/home/ubuntu/AstroCompute_Scripts/initial_clean.py "\
            "/home/ubuntu/data/params.txt /home/ubuntu/"
        objdet_cmd = \
            "/home/ubuntu/miniconda/bin/python " \
            "/home/ubuntu/AstroCompute_Scripts/Aegean_ObjDet.py "\
            "/home/ubuntu/data/params.txt /home/ubuntu/"

        commands = [clone_cmd, clean_cmd, objdet_cmd, timing_cmd]
    else:
        commands = [clone_cmd, timing_cmd]

    mess = {}

    mess['proc_name'] = proc_name
    mess['key_name'] = "data/*"
    mess['bucket'] = proc_name
    mess["command"] = commands

    return json.dumps(mess)


if __name__ == "__main__":
    pass
