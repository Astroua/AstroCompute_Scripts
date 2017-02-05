
from flask_wtf import Form
from wtforms import StringField, BooleanField, PasswordField, IntegerField, \
    TextField, FloatField, FileField, SelectMultipleField
from wtforms.validators import DataRequired, Length, Email, Required

from folder_upload import FolderField
from wtforms import validators


class InputForm(Form):
    email = TextField("Email address",
                      validators=[Required("Please provide a valid email address"),
                                  Length(min=6, message=(u'Email address is too short')),
                                  Email(message=(u'That is not a valid email address.'))])
##
    choice_bool = [("T","T"),("F","F")]
##DATA SET PARAMETERS
    visibility = FolderField('MS Name', validators=[DataRequired()])
    target = StringField('Target Name', validators=[DataRequired()])
    obsDate = StringField("Observation Date", validators=[DataRequired()])
    refFrequency = StringField("Reference Frequency (with units)",
                          validators=[DataRequired()])
    spw_choice = StringField("SPW & Channel selection", default='')
##FLAGS
    choice_run=[("T","T"),("U","U")]
    runClean = SelectMultipleField("Analysis Plane Choice; Only UV Plane (U), only Image Plane or Both (T)",\
        choices = choice_run, default =["T"])
    runObj = SelectMultipleField("Run Source Detection Algorithm?",choices = choice_bool, default =["F"])
    def_times = SelectMultipleField("Manually define start and stop times? (set to F if you want to use full observation)",\
        choices = choice_bool, default =["F"])
##TIMEBINS
    intervalSizeH = IntegerField("Hour Interval Size (allowed values -> 0 to 23)", default=0)
    intervalSizeM = IntegerField("Minute Interval Size (allowed values -> 0 to 59)", default=1)
    intervalSizeS = IntegerField("Second Interval Size (allowed values -> 1e-6 to 59)", default=0)
    startTimeH = StringField("UTC Hour Start Time (specify if def_times=T)", default='')
    startTimeM = StringField("UTC Minute Start Time (specify if def_times=T)", default='')
    startTimeS = StringField("UTC Second Start Time (specify if def_times=T)", default='')
    endTimeH = StringField("UTC Hour End Time (specify if def_times=T)", default='')
    endTimeM = StringField("UTC Minute End Time (specify if def_times=T)", default='')
    endTimeS = StringField("UTC Second End Time (specify if def_times=T)", default='')
##CLEAN PARAMETERS
    imageSize = IntegerField("Image Size", default=256)
    numberIters = IntegerField("# of CLEAN iterations", default=5000)
    cellSize = StringField("Pixel Size (with units)", default='0.3arcsec')
    taylorTerms = IntegerField("# of Taylor Terms", default=1)
    myStokes = StringField("Stokes Axis", default='I')
    thre = StringField("CLEAN Threshold (with units)", default='1mJy')
    outlierFile = FileField("Outlier file for CLEAN")
##OBJECT DETECTION AND SELECTION PARAMETERS-->only need to specify if runObj = 'T'
    seed = IntegerField("Aegean Seed level", default=5)
    flood = IntegerField("Aegean flood level", default=4)
    choice_tels=[("VLA","VLA"),("SMA","SMA"),("NOEMA","NOEMA"),("OTHER","OTHER")]
    tele = SelectMultipleField("Aegean telescope selection",choices = choice_tels, default =["VLA"])
    lat = StringField("Aegean telescope latitude (specify if telescope = OTHER)", default='')
##IMAGE FITTING PARAMETERS-->only need to specify if runClean=T
    choice_mask=[("box","specified targetBox"),("file","CASA .mask file"),("aegean","box around all aegean detected objects")]
    mask_option = SelectMultipleField("CLEAN mask selection",choices = choice_mask, default =["box"])
    mask_file = StringField("CLEAN mask file (specify if mask_option = file)", default='')
    targetBox = StringField("Location of target (pixels; needs to be set if runObj=F)",default='')
    fix_pos = SelectMultipleField("Fix xy position of source for all time-bins?",\
        choices = choice_bool, default =["F"])
    do_monte = SelectMultipleField("Use monte carlo sampling when fitting? (specify if fix_pos=T)",\
        choices = choice_bool, default =["F"])
    choice_integ=[("B","both"),("T","integrated flux"),("F","peak flux")]
    integ_fit = SelectMultipleField("Flux component to record",choices = choice_integ, default =["F"])
##DATA PRODUCT PARAMETERS
    choice_units_f=[("","Jy"),("m","mJy"),("u","uJy")]
    choice_units_t=[("H","hours"),("M","minutes"),("S","seconds")]
    lc_scale_unit = SelectMultipleField("Flux Units",choices = choice_units_f, default =["m"])
    lc_scale_time = SelectMultipleField("Time Units",choices = choice_units_t, default =["M"])
##UV FITTING PARAMETERS
    uv_fit = SelectMultipleField("Fit in the UV plane in addition to Image Plane? (only specify if runClean=T)",\
        choices = choice_bool, default =["F"])
    uv_fix = SelectMultipleField("Fix xy position of source for all time-bins in UV fitting?",\
        choices = choice_bool, default =["F"])
##VARIABILITY ANALYSIS
    var_anal = SelectMultipleField("Do variability analysis?",choices = choice_bool, default =["T"])
    power_spec = SelectMultipleField("Make a power spectrum?",choices = choice_bool, default =["T"])



class LoginForm(Form):
    inputid = StringField('Input ID', validators=[DataRequired()])
    passwd = PasswordField("Password", validators=[DataRequired()])



class EnterDBInfo(Form):
    username = TextField(label='enter username', description="db_enter", validators=[validators.required(), validators.Length(min=0, max=128, message=u'Enter 128 characters or less')])    
    passwd = TextField(label='enter password', description="db_enter_pw", validators=[validators.required(), validators.Length(min=0, max=128, message=u'Enter 128 characters or less')])    
    institute = TextField(label='enter insitution', description="db_enter_in", validators=[validators.required(), validators.Length(min=0, max=128, message=u'Enter 128 characters or less')])
    email_ad = TextField(label='enter email', description="db_enter_email", validators=[validators.required(), validators.Length(min=0, max=128, message=u'Enter 128 characters or less')])
class RetrieveDBInfo(Form):
    username_ret = TextField(label='username to retrieve', description="db_get", validators=[validators.required(), validators.Length(min=0, max=128, message=u'Enter 128 characters or less')])
    passwd_ret = TextField(label='password to retrieve', description="db_get_pw", validators=[validators.required(), validators.Length(min=0, max=128, message=u'Enter 128 characters or less')])



