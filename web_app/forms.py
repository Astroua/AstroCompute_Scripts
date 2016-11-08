
from flask_wtf import Form
from wtforms import StringField, BooleanField, PasswordField, IntegerField, \
    TextField, FloatField, FileField
from wtforms.validators import DataRequired, Length, Email, Required

from folder_upload import FolderField


class InputForm(Form):
    email = TextField("Email address",
                      validators=[Required("Please provide a valid email address"),
                                  Length(min=6, message=(u'Email address is too short')),
                                  Email(message=(u'That is not a valid email address.'))])
    filename = FolderField('File Name', validators=[DataRequired()])
    source_detect = BooleanField('Enable source detection', default=False)
    runObj = source_detect
    target_name = StringField('Name of target', validators=[DataRequired()])
    obsdate = StringField("Observation Date", validators=[DataRequired()])
    reffreq = StringField("Reference Frequency (with units!)",
                          validators=[DataRequired()])

    # visibility = 'v404_jun22_B_Cc7_bp.ms'
    # visibility_uv='v404_jun22_B_Cc7_bp.ms'
    imageSize = IntegerField("Image size", default=6000)
    numberIters = IntegerField("Clean iterations", default=5000)
    cellSize = StringField("Cell size", default='0.02arcsec')
    taylorTerms = IntegerField("Clean Taylor terms", default=1)
    myStokes = StringField("Stokes terms", default='I')
    thre = StringField("Clean threshold", default='10mJy')
    spw_choice = StringField("Clean SPW & Channel select", default='0~7:5~58')
    # only for Aegean_ObjDet.py
    seed = IntegerField("Aegean Seed", default=134656956)
    flood = IntegerField("Aegean flood level", default=4)
    # Define custom requirement.
    tele = StringField("Telescope", default='VLA')
    # only set if not on list in scope2lat func of aegean.py
    lat = StringField("Latitude", default='')  # if no telescope given
    # only for casa_timing_script.py
    ind = IntegerField("Nth detection to run timing on",
                       default=1)

    # only set if not doing OD
    targetBox = StringField("Box around object",
                            default='2982,2937,2997,2947')

    # Annulus size to define source vs. bkg properties
    annulus_rad_inner = FloatField("Inner radius", default=10.)
    annulus_rad_outer = FloatField("Outer radius", default=20.)

    # Run timing script on cutout of image.
    cutout = BooleanField("Cutout region around object", default=True)
    pix_shift_cutout = IntegerField("Pix shift???", default=20)

    # Define an internal parameter for "big data"
    big_data = 'T'

    outlierFile = FileField("Clean outlier file")

    # Why does this need to be set? Or can you also feed it an image?
    runClean = 'T'

    fit_cutout = 'F'
    fix_pos = 'T'

    # Monte Carlo settings
    do_monte = BooleanField("Monte Carlo analysis", default=False)
    nsim = IntegerField("Monte Carlo iterations", default=100)

    # Unsure
    par_fix = 'xyabp'
    integ_fit = 'B'
    uv_fit = 'F'
    do_monte_uv = 'F'
    uv_fix = 'F'

    # for uv fit
    stokes_param = StringField("UV multifit Stokes param", default='I')

    # Set time periods to fit.
    def_times = BooleanField("Define start and stop times to fit",
                             default=False)

    # Interval sizes
    intervalSizeH = IntegerField("Hour Interval Size", default=0)
    intervalSizeM = IntegerField("Minute Interval Size", default=0)
    intervalSizeS = IntegerField("Second Interval Size", default=2)

    # only set if def_times='T'
    startTimeH = StringField("Hour Start Time", default='')
    startTimeM = StringField("Minute Start Time", default='')
    startTimeS = StringField("Second Start Time", default='')
    endTimeH = StringField("Hour End Time", default='')
    endTimeM = StringField("Minute End Time", default='')
    endTimeS = StringField("Second End Time", default='')

    # This doesn't appear to be used in the timing scripts.
    rem_int = 5

    # Time units
    lc_scale_unit = 'm'
    lc_scale_time = 'M'

    opt_clean = StringField("Optimize clean parameters", default='F')


class LoginForm(Form):
    inputid = StringField('Input ID', validators=[DataRequired()])
    passwd = PasswordField("Password", validators=[DataRequired()])


class RegisterForm(Form):
    username = StringField('User name', validators=[DataRequired()])
    email = \
        TextField("Email address",
                  validators=[Required("Please provide a valid email"
                                       "address"),
                              Length(min=6, message=(u'Email address is too'
                                                     ' short')),
                              Email(message=(u'That is not a valid email '
                                             'address.'))])
    institution = StringField("Institution", validators=[DataRequired()])
    passwd = PasswordField("Password", validators=[DataRequired()])
