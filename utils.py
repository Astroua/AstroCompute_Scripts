
import json
import re


def convert_param_format(filename, to="json"):
    '''
    Convert the text format into a dictionary, or JSON.
    '''

    with open(filename, "r") as f:
        params = f.read()

    params = params.split("\n")
    param_dict = {}

    for param in params:
        # Skip comment lines
        if len(param) == 0:
            continue

        if param[0] == "#":
            continue

        # Strip whitespace
        param = re.sub("\s", "", param)

        splits = param.split("=")

        key = splits[0]
        value = splits[1]

        if len(value.split("#")) > 1:
            value = value.split("#")[0]

        param_dict[key] = value.strip("'")

    if to is "json":
        return json.dumps(param_dict)
    else:
        return param_dict


def load_json(filename):
    '''
    Load in a JSON formatted text file.
    '''

    with open(filename, "r") as f:
        contents = json.load(f)

    for key in contents:
        contents[key] = str(contents[key])

    return contents


def is_power2(num):
    ''' Check if input imsize is a power of 2^n, in order to
    optimize cleaning speed

    num: image size in pixels

    return: True/False
    '''
    return(num != 0 and ((num & (num-1)) == 0))
