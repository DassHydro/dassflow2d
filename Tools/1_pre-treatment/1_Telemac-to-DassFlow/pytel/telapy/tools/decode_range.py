"""
   Function to handle range strings
"""
from utils.exceptions import TelemacException


def decode_range(string):
    """
    Transform a string in format [1,2:8,..,4] into a list
    """
    res = []
    # Checking that beginning and end of the string are [ ]
    if string[0] != '[' and string[:-1] != "]":
        raise TelemacException("Invalid range format for %s" % string)

    # Splitting values that should be separated by commas
    tmp_list = string[1:-1].split(",")

    for item in tmp_list:
        # Wide range item
        if ":" in item:
            i, j = item.split(":")
            for val in range(int(i), int(j)+1):
                res.append(val)
        # Isolated item
        else:
            res.append(int(item))

    return res
