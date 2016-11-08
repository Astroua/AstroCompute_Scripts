

from flask import Flask, request, redirect, url_for, flash, \
    render_template, send_from_directory, jsonify, request
from werkzeug.utils import secure_filename
from astropy import log
from flask_bootstrap import Bootstrap
from flask_login import login_required, current_user

from aws_controller.upload_download_s3 import upload_to_s3

from web_app import InputForm, LoginForm, RegisterForm, ResultInfo


UPLOAD_FOLDER = 'uploads/'

app = Flask(__name__)
Bootstrap(app)
app.config.from_object('config.AWSConfig')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

key = app.config['AWS_KEY']
secret = app.config['AWS_SECRET']


@app.route("/")
def default():
    # display welcome page
    return render_template('index.html')


@app.route('/submit', methods=['GET', 'POST'])
@login_required
def submit():
    form = InputForm()
    form.filename(multiple="")
    if request.method == "POST":
        if form.validate_on_submit():
            return render_template('job_submission.html')
        else:
            return render_template('submit.html', form=form)
    return render_template('submit.html', form=form)


@app.route('/upload', methods=['POST'])
def upload(timestamp):
    if request.method == 'POST':
        data_file = request.files('file')
        file_name = secure_filename(data_file.filename)
        upload_to_s3(timestamp, file_name, create_bucket=True,
                     aws_access={"aws_access_key_id": key,
                                 "aws_secret_access_key": secret})

        return jsonify(name=file_name)


@app.route("/login", methods=['GET', 'POST'])
def login():
    '''
    Login form for users
    '''
    form = LoginForm()
    if form.validate_on_submit():
        flash("Logging in: " + form.inputid.data)
        return redirect("/")
    return render_template('login.html', title='Login', form=form)


@app.route("/register", methods=['GET', 'POST'])
def register():
    '''
    Register form for new users
    '''
    form = LoginForm()
    if form.validate_on_submit():
        flash("Logging in: " + form.inputid.data)
        return redirect("/")
    return render_template('register.html', title='Register', form=form)


@app.route("/summary", methods=['POST'])
@login_required
def summary(params):
    output = ResultInfo(params)
    job_name = params["job_name"]

    return render_template('summary.html',
                           title='Results for {}'.format(job_name),
                           output=output)


if __name__ == '__main__':
    log.setLevel(10)
    app.run(debug=True)
