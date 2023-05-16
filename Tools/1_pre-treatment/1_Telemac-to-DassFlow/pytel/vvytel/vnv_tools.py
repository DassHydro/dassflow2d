
"""
Functions for validation that compare resultat and compute norms
"""
# _____                               ______________________________
# ____/ Test case validation Toolkit /_____________________________/
#
#
import numpy as np
from data_manip.extraction.telemac_file import TelemacFile
from utils.exceptions import TelemacException

def norm_l1(diff, mass=None):
    """
    Compute the l1 norm for the difference between v1 and v2

    @param diff (np.array) Vector containing (x1 - x2)
    @param mass (np.array) Mass factor (default: None)
    """
    if mass is None:
        return np.sum(abs(diff))/len(diff)
    else:
        return np.sum(mass*abs(diff))/np.sum(mass)

def norm_l2(diff, mass=None):
    """
    Compute the l2 norm for the difference between v1 and v2

    @param diff (np.array) Vector containing (x1 - x2)
    @param mass (np.array) Mass factor (default: None)
    """
    if mass is None:
        return np.sqrt(np.sum(abs(diff*diff))/len(diff))
    else:
        return np.sqrt(np.sum(mass*abs(diff*diff))/np.sum(mass))

def norm_linf(diff):
    """
    Compute the linf norm for the difference between v1 and v2

    @param diff (array) Vector containing (x1 - x2)
    """
    return np.max(abs(diff))

def compute_norm(diff, norm='linf', mass=None):
    """
    Compute norm of a vector

    @param diff (array) Vector containing (x1 - x2)
    @param mass (array) Mass factor (default: None)
    @param norm (str) type of norm
    """
    if norm == 'l1':
        return norm_l1(diff, mass=mass)
    elif norm == 'l2':
        return norm_l2(diff, mass=mass)
    else:
        return norm_linf(diff)

def compute_diff(data1, data2, relative=False):
    """
    Compute difference of two vectors

    @param data1 (array) usually the model or computation
    @param data2 (array) usually the reference data
    @param relative (bool) use relative difference
    """
    if relative:
        diff = data1 - data2
        for i, _ in enumerate(diff):
            if abs(diff[i]) > 1e-42:
                a_max = max(abs(data1[i]), abs(data2[i]))
                diff[i] = diff[i]/a_max
            else:
                diff[i] = 0.0
    else:
        diff = data1 - data2
    return diff

def check_compatibility(res1, res2, record, check_name=True):
    """
    Check if two TelemacFile object contain the same variables,
    record times and number of points

    @param res1 (TelemacFile) first telemac-mascaret result file
    @param res2 (TelemacFile) second telemac-mascaret result file
    @param record (int) Record to check
    @param check_name (bool) If False only check that we have the same number
    of variable but not necessaraly the same name
    """
    # Checking that we have the same variables
    if check_name:
        if res1.varnames != res2.varnames:
            raise TelemacException(\
                    "Not the same variables in each files\n"
                    "Variables in {}: {}\nVariables in {}: {}"\
                    .format(res1.file_name, res1.varnames,
                            res2.file_name, res2.varnames))
    else:
        if len(res1.varnames) != len(res2.varnames):
            raise TelemacException(\
                "Not the same variables in each files\n"
                "Number of variables in {}: {}\nNumber of variables in {}: {}"\
                .format(res1.file_name, len(res1.varnames),
                        res2.file_name, len(res2.varnames)))

    # Checking that we have the same record time
    time1 = res1.get_data_time(record)
    time2 = res2.get_data_time(record)

    if abs(time1-time2) > 1e-6:
        raise TelemacException(\
                "Not the same time in each file\n"
                "In {}: {}\nIn {}: {}".format(res1.file_name, time1,
                                              res2.file_name, time2))

    # Checking that we have the same number of points
    if res1.npoin3 != res2.npoin3:
        raise TelemacException(\
                "Not the same time in each file\n"
                "In {}: {}\nIn {}: {}"\
                .format(res1.file_name, res1.npoin3,
                        res2.file_name, res2.npoin3))
