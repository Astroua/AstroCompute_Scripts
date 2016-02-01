

from flask import Flask, request, redirect, url_for, \
    render_template, send_from_directory, jsonify
from astropy import log
from flask_wtf import Form
from wtforms import StringField, BooleanField
from wtforms.validators import DataRequired

UPLOAD_FOLDER = 'uploads/'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'cheese'

key = app.config['S3_KEY']
secret = app.config['S3_SECRET']
bucket = app.config['S3_BUCKET']


class MyForm(Form):
    name = StringField('name', validators=[DataRequired()])


@app.route("/")
def default():
    # display welcome page
    return render_template('index.html')


@app.route('/submit', methods=('GET', 'POST'))
def submit():
    form = MyForm()
    if form.validate_on_submit():
        return redirect('/success')
    return render_template('submit.html', form=form)


if __name__ == '__main__':
    log.setLevel(10)
    app.run(debug=True)
