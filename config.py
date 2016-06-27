
from boto import Config
import os


class AWSConfig(object):
    """docstring for AWSConfig"""
    DEBUG = True

    config = Config()
    config.load_credential_file(os.path.join(os.path.expanduser("~"),
                                             ".aws/credentials"))
    info = config.items("default")[2:]
    AWS_KEY = info[0][1]
    AWS_SECRET = info[1][1]

    WTF_CSRF_ENABLED = True
    SECRET_KEY = 'cheese'
