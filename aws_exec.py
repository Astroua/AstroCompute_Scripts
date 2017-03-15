
import json
from time import sleep, time

from boto import sqs
from boto import ses

from aws_controller.controller import WORKER_SCRIPT
from aws_controller.launch_instance import launch
from aws_controller.upload_download_s3 import set_bucket_lifetime


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

    proc_name = params['target'].lower() + "_" + timestamp
    key = aws_settings['key']
    secret = aws_settings['secret']
    region = aws_settings['region']

    # Create the queue and submit the message.
    queue_name = proc_name
    resp_queue_name = "resp_" + proc_name

    queue = sqs.connect_to_region(region,
                                  aws_access_key_id=key,
                                  aws_secret_access_key=secret).create_queue(queue_name)

    msg = json_message(params, proc_name)
    queue.write(queue.new_message(body=msg))

    # Launch instance

    user_data = WORKER_SCRIPT % \
        {"USER": 'ubuntu',
         "QUEUE_NAME": proc_name,
         "KEY": key,
         "SECRET": secret,
         "REGION": region,
         "RESP_QUEUE_NAME": resp_queue_name,
         "CUSTOM_LINES": "/usr/bin/git clone https://github.com/Astroua/AstroCompute_Scripts.git"}
         # "CUSTOM_LINES": 'ssh-keyscan github.com >> /home/ubuntu/.ssh/known_hosts\nsu - ubuntu -c "/usr/bin/git clone git@github.com:e-koch/AstroCompute_Scripts.git"\nexport USER="ubuntu"'}

    inst = launch(key_name=None, region=region,
                  image_id=aws_settings['image_id'],
                  instance_type=aws_settings['instance_type'],
                  security_groups=aws_settings['security_groups'],
                  initial_check=False, user_data=user_data)

    # Wait for return message
    # Total wait time set to 1 day
    time_max = 3600 * 24
    time_wait = 200
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
        with open("test_response.txt", "w") as f:
            json.dump(content, f)

    # Clean-up queues
    queue.delete()
    resp_queue.delete()

    # Send out email to the user.
    email_service = ses.connect_to_region(region,
                                          aws_access_key_id=key,
                                          aws_secret_access_key=secret)

    # The parameters should have the user's database values.
    email_address = params["user_email"]
    username = params["user_name"]
    time_limit = params["time_limit"]
    summary_url = params["summary_url"]
    project_email = params["project_email"]

    subject = "AstroCompute Timing Script Job: {} Completion".format(proc_name)
    body = email_body(username, proc_name, summary_url, time_limit,
                      project_email)
    email_service.send_email(project_email, subject, body, email_address,
                             format='text')

    # Set the expiration on the S3 bucket.
    set_bucket_lifetime(proc_name, days=time_limit,
                        aws_access={"aws_access_key_id": key,
                                    "aws_secret_access_key": secret})

    # Return success/failure

    # Add in the job name to parameters
    params["job_name"] = proc_name

    return params


def json_message(params, proc_name):
    '''
    Create the json message to feed to the worker instance.
    '''

    # Commands to run
    chmod_cmd = \
        "chmod -R 777 /home/ubuntu/data/" + params['visibility']

    timing_cmd = \
        "/home/ubuntu/casa-release-4.7.1-el7/bin/casa --nologger " \
        "--logfile data_products/" + proc_name + "_timing.log -c "\
        "/home/ubuntu/AstroCompute_Scripts/casa_timing_script.py "\
        "/home/ubuntu/data/params.txt /home/ubuntu/ "\
        "/home/ubuntu/AstroCompute_Scripts/"

    # Check whether Object detection should be run
    if params['runObj'] == "T":
        clean_cmd = \
            "/home/ubuntu/casa-release-4.7.1-el7/bin/casa --nologger " \
            "--logfile data_products/" + proc_name + "_initclean.log -c "\
            "/home/ubuntu/AstroCompute_Scripts/initial_clean.py "\
            "/home/ubuntu/data/params.txt /home/ubuntu/"
        objdet_cmd = \
            "/home/ubuntu/miniconda/bin/python " \
            "/home/ubuntu/AstroCompute_Scripts/Aegean_ObjDet.py "\
            "/home/ubuntu/data/params.txt /home/ubuntu/"

        # commands = [clean_cmd, objdet_cmd, timing_cmd]
        commands = [chmod_cmd, clean_cmd, objdet_cmd, timing_cmd]
    else:
        commands = [chmod_cmd, timing_cmd]

    params['proc_name'] = proc_name
    params['key_name'] = "data/*"
    params['bucket'] = proc_name
    params["command"] = commands

    return json.dumps(params)


def email_body(username, job_name, url, time_limit, project_email):
    '''
    Generate string for email body.
    '''

    email_str = \
        '''
        Hello USER,\n
        \n
        Your timing script job JOB is complete. The analysis results may be
        downloaded from the link below. Note that data products are
        automatically deleted after TIME days. Please retrieve the data
        products as soon as possible!\n
        \n
        URL \n
        \n
        Questions or issues should be sent to EMAIL_HELP.\n
        '''

    email_str = email_str.replace("USER", username)
    email_str = email_str.replace("JOB", job_name)
    email_str = email_str.replace("URL", url)
    email_str = email_str.replace("TIME", str(time_limit))
    email_str = email_str.replace("EMAIL_HELP", project_email)

    return email_str


if __name__ == "__main__":
    pass
