
import imp


def getVar(filename):
    '''
    Easy way to read in parameters from file
    '''
    with open(filename) as f:
        data_params = imp.load_source('data_params', '', f)

    return data_params


def is_power2(num):
    ''' Check if input imsize is a power of 2^n, in order to
    optimize cleaning speed

    num: image size in pixels

    return: True/False
    '''
    return(num != 0 and ((num & (num-1)) == 0))
