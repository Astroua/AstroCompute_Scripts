
from flask_wtf import Form
from wtforms import StringField, BooleanField
from wtforms.validators import DataRequired


class InputForm(Form):
    filename = StringField('File Name', validators=[DataRequired()])
    source_detect = BooleanField('Enable source detection', default=False)


class LoginForm(Form):
    inputid = StringField('Input ID', validators=[DataRequired()])
    passwd = StringField("Password", validators=[DataRequired()])
