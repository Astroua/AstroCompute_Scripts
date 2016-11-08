
from flask_wtf import Form
from wtforms import StringField, BooleanField, PasswordField, IntegerField, \
    TextField, FloatField, FileField, SelectMultipleField
from wtforms.validators import DataRequired, Length, Email, Required

from folder_upload import FolderField



class InputForm(Form):
    email = TextField("Email address",
                      validators=[Required("Please provide a valid email address"),
                                  Length(min=6, message=(u'Email address is too short')),
                                  Email(message=(u'That is not a valid email address.'))])
    filename = FolderField('File Name', validators=[DataRequired()])
    source_detect = BooleanField('', default=False)
    runObj = source_detect
    target_name = StringField('Name of target', validators=[DataRequired()])
    obsdate = StringField("Observation Date", validators=[DataRequired()])
    reffreq = StringField("Reference Frequency (with units!)",
                          validators=[DataRequired()])
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
    choice_tels=[("VLA","VLA"),("SMA","SMA"),("NOEMA","NOEMA"),("OTHER","OTHER")]
    tele = SelectMultipleField("Telescope",choices = choice_tels, default =["VLA"])
    # only set if not on list in scope2lat func of aegean.py
    lat = StringField("Latitude (specify if telescope = OTHER)", default='')  # if no telescope given
    # only for casa_timing_script.py
    #ind = IntegerField("Nth detection to run timing on",
    #                   default=1) **enable multiple objects later**

    # only set if not doing OD
    mask_option = BooleanField('', default=False)
    targetBox = StringField("Size of box around object (pixels)",
                            default='2982,2937,2997,2947')

    # Annulus size to define source vs. bkg properties
    annulus_rad_inner = FloatField("Inner annulus radius (pixels)", default=10.)
    annulus_rad_outer = FloatField("Outer annulus radius (pixels)", default=20.)

    # Run timing script on cutout of image.
    cutout = BooleanField("Cutout region around object", default=True)
    pix_shift_cutout = IntegerField("Size of cutout image (pixels from target box)", default=20)

    # Define an internal parameter for "big data"
    big_data = BooleanField('delete original image to save space?', default=True)

    outlierFile = FileField("Clean outlier file")

    # Why does this need to be set? Or can you also feed it an image?
    runClean = BooleanField('run clean?', default=True)

    #why is this here?
    fit_cutout = BooleanField('re-fit each time bin image with no cleaning?', default=False)

    fix_pos = BooleanField('fix position in individual images?', default=True)

    # Monte Carlo settings
    do_monte = BooleanField("Use monte carlo sample position of full data set", default=False)
    nsim = IntegerField("Monte Carlo iterations", default=100)

    # Unsure
    #par_fix = 'xyabp' wtf alex??
    #integ_fit = 'B' wtf alex?? if you want a lsit of options use SelectMultipleField?
    uv_fit = BooleanField("uv plane fitting??", default=False)
    do_monte_uv = BooleanField("monte carlo sample fixed position in uv fitting?", default=False)
    uv_fix = BooleanField("fix position in uv fitting?", default=False)

    # for uv fit
    stokes_param = StringField("stokes parameter for uv fitting", default='I')

    # Set time periods to fit.
    def_times = BooleanField("Manually define start and stop times to fit to:",
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

    # check amount of remaining time after dividing time bins
    rem_int = IntegerField("if remaining time >, include in time intervals", default=5)

    # Time units
    choice_units_f=[("","Jy"),("m","mJy"),("u","uJy")]
    choice_units_t=[("H","hours"),("M","minutes"),("S","seconds")]
    lc_scale_unit = SelectMultipleField("Flux Units",choices = choice_units_f, default =["m"])
    lc_scale_time = SelectMultipleField("Time Units",choices = choice_units_t, default =["M"])

    opt_clean = StringField("Optimize clean parameters", default='F') #wtf alex??


class LoginForm(Form):
    inputid = StringField('Input ID', validators=[DataRequired()])
    passwd = PasswordField("Password", validators=[DataRequired()])
