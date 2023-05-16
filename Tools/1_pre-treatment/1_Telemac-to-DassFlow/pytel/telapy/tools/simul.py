# -*- coding: utf-8 -*-
"""
Definition of the cost function for 'genop' optimizer

Author(s) : Fabrice Zaoui

Copyright EDF 2016-2017

:param 'x': array of optimization variables
:return: value of the cost function
"""
import numpy as np


# user function
def sim_telemac(varx):
    """
    The cost function to be written by user
    """
    val = varx
    return val


# test function
def rosenbrock(varx):
    """
    https://en.wikipedia.org/wiki/Rosenbrock_function
    """
    val = (1.-varx[0])**2+100.*(varx[1]-varx[0]**2)**2
    return val


# test function
def maccormick(varx):
    """
    https://en.wikipedia.org/wiki/Test_functions_for_optimization
    """
    val = np.sin(varx[0]+varx[1])+(varx[0]-varx[1])**2 - \
        1.5*varx[0]+2.5*varx[1]+1.
    return val


# test function
def easom(varx):
    """
    https://en.wikipedia.org/wiki/Test_functions_for_optimization
    """
    val = -np.cos(varx[0])*np.cos(varx[1]) * \
        np.exp(-((varx[0]-np.pi)**2+(varx[1]-np.pi)**2))
    return val


# test function
def eggholder(varx):
    """
    https://en.wikipedia.org/wiki/Test_functions_for_optimization
    """
    val = -(varx[1]+47.)*np.sin(np.sqrt(np.abs(0.5*varx[0]+varx[1]+47.))) - \
        varx[0]*np.sin(np.sqrt(np.abs(varx[0]-varx[1]-47.)))
    return val


# test function
def linear(varx):
    """
    a simple linear function of three variables
    """
    val = varx[0] + varx[1] - varx[2]
    return val


# test function
def abslinear(varx):
    """
    a simple absolute value function
    """
    val = np.abs(varx[0]-2.)
    return val
