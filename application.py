

from flask import Flask, request, redirect, url_for, \
    render_template, send_from_directory, jsonify, request
from astropy import log
from flask_wtf import Form
from wtforms import StringField, BooleanField
from wtforms.validators import DataRequired

from aws_controller.upload_download_s3 import upload_to_s3


UPLOAD_FOLDER = 'uploads/'

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'cheese'

key = app.config['S3_KEY']
secret = app.config['S3_SECRET']


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


@app.route('/upload', methods=['POST'])
def upload(timestamp):
    if request.method == 'POST':
        data_file = request.files.get('file')
        file_name = data_file.filename
        upload_to_s3(timestamp, file_name, create_bucket=True,
                     aws_access={"aws_access_key_id": key,
                                 "aws_secret_access_key": secret})

        # return jsonify(name=file_name)
        return jsonify(name=file_name)

if __name__ == '__main__':
    log.setLevel(10)
    app.run(debug=True)
