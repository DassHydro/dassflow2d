"""
Functions for validation that compare resultat and compute norms
"""
# _____                                         ______________________________
# ____/ Test case validation Toolkit /_____________________________/
#
#
from numpy import sqrt, zeros
from utils.exceptions import TelemacException
###############################################################################
def mapdiff(a_1, a_2, notime=False, noname=False, relative=False):
    """
    Function used in the xml files to do the comparaison between two files.
    @brief Create a new values containing the diff of the values
    of the two arguments

    @param a_1 A Values object containing the keys {support,names,values,time}
    @param a_2 A Values object containing the keys {support,names,values,time}
    @param notime (optional) allow comparison of result with different time
    @param noname (optional) Not checking that names are the same
    @param relative (boolean) If true relative error is computed otherwise
                              absolute one


    @return A 4-uple (time,name_of_the_variables,support,values)
    """
    # Cheking they have the same shape should be (ntime,nvar,npoin)
    diff = zeros(a_2.values.shape, dtype=a_2.values.dtype)
    # With ntime = 1
    if a_1.values.shape != a_2.values.shape:
        raise TelemacException(\
              "Error in files the two array do not have the same shape"\
              " (ntimestep, nvar, npoin)\n"+
              "a_1 shape: " + repr(a_1.values.shape) + "\n"
              "a_2 shape: " + repr(a_2.values.shape))
    # Shape of the values should be (ntime,nvar,npoin)
    _, nvar, npoin = a_2.values.shape
    # Checking that it is the same variables in each files unless nonam is true
    if not noname:
        for ivar in range(nvar):
            if a_1.names[ivar][0:15] != a_2.names[ivar][0:15]:
                raise TelemacException(\
                    "Could not found '{}' from a_1 in a_2:\n{}".format( \
                    a_1.names[ivar][0:16], a_2.names))
    # Checking if we have the same time step
    if not notime:
        if abs(a_1.time[0] - a_2.time[0]) > 1e-6:
            raise TelemacException(\
                'The time of the two times are different \n' + \
                str(a_1.time[0]) + ' for a_1\n' + \
                str(a_2.time[0]) + ' for a_2')
    # Making a_1-a_2 for each variable and each point
    # TODO: Optimize using numpy computation (remove loops)
    if relative:
        for ivar in range(nvar):
            for i in range(npoin):
                ldiff = a_2.values[0][ivar][i] \
                           - a_1.values[0][ivar][i]
                if abs(ldiff) > 1e-42:
                    a_max = max(abs(a_2.values[0][ivar][i]), \
                                abs(a_1.values[0][ivar][i]))
                    diff[0][ivar][i] = ldiff/a_max
                else:
                    diff[0][ivar][i] = 0.0
    else:
        for ivar in range(nvar):
            for i in range(npoin):
                diff[0][ivar][i] = a_2.values[0][ivar][i] \
                                  - a_1.values[0][ivar][i]

    return a_2.time, a_2.names, a_2.support, diff
###############################################################################
def norm_l1(diff):
    """
    Compute the l1 norm for the difference between v1 and v2

    @param diff Vector containing (x1 - x2)
    """
    err = sum(abs(diff))/len(diff)
    return err
###############################################################################
def norm_l2(diff):
    """
    Compute the l2 norm for the difference between v1 and v2

    @param diff Vector containing (x1 - x2)
    """
    err = sqrt(sum(abs(diff*diff)))/len(diff)
    return err
###############################################################################
def norm_linf(diff):
    """
    Compute the linf norm for the difference between v1 and v2

    @param diff Vector containing (x1 - x2)
    """
    err = max(abs(diff))
    return err
###############################################################################
def checkval(a_0, eps, norm='linf'):
    """
    @brief Will loop on all variable and display the max error

    @param a_0 A Values object containg the difference between two results
    @param eps The epsilon for each variable or a global one

    @param norm (str) Norm to apply (linf, l1, l2)
    @return True if all the variable max are below EPS; False otherwise
    """
    _, nvar, _ = a_0.values.shape
    # Getting eps for each variable
    force_success = False
    if not eps:
        force_success = True
        # if eps is empty setting default value i.e. 1.e-4
        print(" "*8+"~> TODO: Set Epsilon value")
        final_eps = [1.e-4]*nvar
    elif len(eps) == 1:
        # If only one was given using it for all the variables
        final_eps = [eps[0]] * nvar
    elif len(eps) != nvar:
        # Otherwise one mus be given for each variable
        raise TelemacException(\
            "Error in length of espilon is {} should be {}"\
            .format(len(eps), nvar))
    else:
        final_eps = eps
    print(" " * 8 + "~> Validation for time (in reference file):" + \
            str(a_0.time[0]))
    # Loop on variables
    value = False
    for ivar in range(nvar):
        if norm == 'l1':
            err = norm_l1(a_0.values[0][ivar])
        elif norm == 'l2':
            err = norm_l2(a_0.values[0][ivar])
        elif norm == 'linf':
            err = norm_linf(a_0.values[0][ivar])
        else:
            err = norm_linf(a_0.values[0][ivar])
        print(" " * 10 + "- Difference for variable " + a_0.names[ivar] + \
                ": " + str(err))
        if err >= final_eps[ivar]:
            print(" " * 12 + "Epsilon reached" + str(final_eps[ivar]))
            value = True
    # Printing empty line
    print("")
    return value and (not force_success)
