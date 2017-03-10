from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
from forms import InputForm
from flask_bootstrap import Bootstrap
from flask_login import LoginManager




UPLOAD_FOLDER = 'uploads/'

app = Flask(__name__)
Bootstrap(app)
app.config.from_object('config.AWSConfig')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
db = SQLAlchemy(app)

key = app.config['AWS_KEY']
secret = app.config['AWS_SECRET']

from app import models
from results_class import ResultInfo
