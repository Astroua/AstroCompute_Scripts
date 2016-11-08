
from boto import Config
import os


class AWSConfig(object):
    """docstring for AWSConfig"""
    DEBUG = True

    # edit the URI below to add your RDS password and your AWS URL
    # The other elements are the same as used in the tutorial
    # format: (user):(password)@(db_identifier).amazonaws.com:3306/(db_name)

    SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://ualberta:H1743m322@user-database.cpluopkax7lm.us-west-2.rds.amazonaws.com:3306/user_db'
    SQLALCHEMY_POOL_RECYCLE = 3600

    config = Config()
    config.load_credential_file(os.path.join(os.path.expanduser("~"),
                                             ".aws/credentials"))
    info = config.items("default")[2:]
    AWS_KEY = info[0][1]
    AWS_SECRET = info[1][1]

    WTF_CSRF_ENABLED = True
    SECRET_KEY = 'cheese'
