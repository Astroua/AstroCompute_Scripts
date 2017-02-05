
import os 
from utils import is_power2
from astropy import units as u

def test_target(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_obsDate(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_refFrequency(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_intervalSizeH(form,field):
	if field.data not in xrange(0,24):
		raise ValueError("Value must be between 0 and 23. Please change input.")
def test_intervalSizeM(form,field):
	if field.data not in xrange(0,60):
		raise ValueError("Value must be between 0 and 59. Please change input.")
def test_intervalSizeS(form,field):
	if field.data < 0. or field.data >= 60.:
		raise ValueError("Value must be between 0.0 and 59.9. Please change input.")
def test_visibility(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if not os.path.isdir(field.data):
		raise ValueError("This data set does not appear to be in the correct format of a CASA MS. Please check input.")
def test_imageSize(form,field):
	if not isinstance(field.data,int):
		raise ValueError("Value entered is not a integer. Please change input to an integer.")
	if not is_power2(field.data):
		raise ValueError("Value entered is not a power of 2^n. To optimize cleaning process please change input.")
def test_numberIters(form,field):
	if not isinstance(field.data,int):
		raise ValueError("Value entered is not a integer. Please change input to an integer.")
	if field.data > 20000:
		raise ValueError("Value is greater then the maximum allowed niters. Please enter a value < 20000.")
def test_cellSize(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	u.Unit('arcsec')
	letters=re.findall("[a-zA-Z]+",field.data)
	number=re.findall(r"[-+]?\d*\.\d+|\d+", field.data)
	if not u.Unit(letters[0]).is_equivalent:
		raise ValueError("Incorrect unit specified. Please use arcsec.")
def test_taylorTerms(form,field):
	if not isinstance(field.data,int):
		raise ValueError("Value entered is not a integer. Please change input to an integer.")
	if field.data < 1:
		raise ValueError("Value must be 1 or greater. Please change input.")
def test_myStokes(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_thre(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	u.Unit('mJy')
	u.Unit('uJy')
	letters=re.findall("[a-zA-Z]+",field.data)
	number=re.findall(r"[-+]?\d*\.\d+|\d+", field.data)
	if not u.Unit(letters[0]).is_equivalent:
		raise ValueError("Incorrect unit specified. Please use mJy or uJy.")
def test_spw_choice(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	string=field.data
	str_split=string.split(':')
	if '~' not in str_split[0] or '~' not in str_split[1]:
		raise ValueError("Value entered is not in the correct format. Please change input to X~X:X~X format.")
def test_mask_option(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	allowed=['box','mask','aegean']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'box', 'mask', or 'aegean'.")
def test_mask_file(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_var_anal(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_power_spec(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_seed(form,field):
	if not isinstance(field.data,int):
		raise ValueError("Value entered is not a integer. Please change input to an integer.")
def test_flood(form,field):
	if not isinstance(field.data,int):
		raise ValueError("Value entered is not a integer. Please change input to an integer.")
def test_tele(form,field):
	allowed=['VLA','SMA','NOEMA']
	if field.data not in allowed:
		raise ValueError("Value entered is not in the default telescope list. Please enter latitude value in degrees in lat parameter instead.")
def test_lat(form,field):
	if not isinstance(field.data,float):
		raise ValueError("Value entered is not a float. Please change input to an float.")
def test_runObj(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_targetBox(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	string=field.data
	str_split=string.split(',')
	if len(str_split) != 4:
		raise ValueError("Value entered is not in the correct format. Please change input to 'X1,Y1,X2,Y2'.")
def test_outlierFile(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_runClean(form,field):
	allowed=['T','F','U']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F' or 'U'.")
def test_fix_pos(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_do_monte(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_par_fix(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	allowed=['x','y','b','a','p']
	string_list=list(field.data)
	if not set(string_list).issubset(set(allowed)):
		raise ValueError("Value entered is invalid. Elements of string may be x (peak x pos), y (peak y pos), f (peak flux), a (major axis),b (minor axis) or, p (position angle).")
def test_integ_fit(form,field):
	allowed=['B','T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' (integrated only),'F' (peak only), or B' (both).")
def test_uv_fit(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_uv_fix(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_stokes_param(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
def test_def_times(form,field):
	allowed=['T','F']
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'T' or 'F'.")
def test_startTimeH(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if int(field.data) not in xrange(0,24):
		raise ValueError("Value entered is invalid. Please change input to a number between 0 and 24.")
def test_startTimeM(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if int(field.data) not in xrange(0,60):
		raise ValueError("Value entered is invalid. Please change input to a number between 0 and 59.")
def test_startTimeS(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if int(field.data) < 0. or int(field.data) >= 60.:
		raise ValueError("Value entered is invalid. Please change input to a number between 0.0 and 59.9")
def test_endTimeH(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if int(field.data) not in xrange(0,24):
		raise ValueError("Value entered is invalid. Please change input to a number between 0 and 24.")
def test_endTimeM(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if int(field.data) not in xrange(0,60):
		raise ValueError("Value entered is invalid. Please change input to a number between 0 and 59.")
def test_endTimeS(form,field):
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if int(field.data) < 0. or int(field.data) >= 60.:
		raise ValueError("Value entered is invalid. Please change input to a number between 0.0 and 59.9")
def test_lc_scale_unit(form,field):
	allowed=['','m','u']
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter '' or 'm', or 'u'.")
def test_lc_scale_time(form,field):
	allowed=['H','M','S']
	if not isinstance(field.data,str):
		raise ValueError("Value entered is not a string. Please change input to a string format.")
	if field.data not in allowed:
		raise ValueError("Value entered is invalid. Please enter 'H' or 'M', or 'S'.")