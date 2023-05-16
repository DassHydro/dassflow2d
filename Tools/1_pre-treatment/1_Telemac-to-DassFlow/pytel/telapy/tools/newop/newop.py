#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
    A quasi-Newton optimizer based on the SciPy 'minimizer' function

    Author(s): Fabrice Zaoui

    Copyright EDF 2017-2018

    Comments : using the L-BFGS-B optimization algorithm
"""

from telapy.tools.newop.validate import validate
from telapy.tools.newop.numval import numv
from scipy.optimize import minimize


class Newop(object):
    """
    The base class for the SciPy Opimization
    """

    def __init__(self, d_x=1.e-6, maxfun=2000, verbose=True):
        """
        Initialize some algorithmic parameters to default values
        :return: a new object from Newop
        """
        # value of the step size for the numerical gradient
        self.d_x = d_x
        self.vdx = None
        # maximum number of function evaluations
        self.maxfun = maxfun
        # print information
        self.verbose = verbose
        # name of the cost function
        self.function = ""
        # number of optimization variables
        self.nvar = 0
        # maximum number of calls to the simulation function
        self.maxfun = maxfun
        # lower and upper bounds of optimization variables
        self.bounds = []
        # intialization is not done
        self.__ready = False

    def initialize(self, func, nvar, bounds, vdx=None):
        """
        Description of the minimization problem
        :param 'f': the name of the python function where the cost function
            is implemented (type: str)
        :param 'n': the number of optimization variables (type: int)
        :param 'bounds': a numpy array for the lower and upper bounds of
            optimization variables (type: int)
        :param 'vdx': the finite difference steps for each variable
        :return: error code = 0 if successful (type: int)
        """
        # Test the arguments
        error = validate(func, nvar, bounds, vdx, self.verbose)
        if error:
            if self.verbose:
                print('--> newop Error!')
                print('\tIncorrect arguments.')
            return -1
        # Save all the information
        self.function = func
        self.nvar = nvar
        self.bounds = bounds
        self.vdx = vdx
        # Job done
        self.__ready = True
        return 0

    def optimize(self, varx0, niter=100, nproc=1):
        """
        Do the l-bfgs-b algorithm with 'niter' iterations
        :param 'varx0': initial guess
        :param 'niter': the maximal number of iterations
        :param 'nproc': the number of processes (parallel simulations)
        :return: a tuple of two elements at the end of the optim process:
                    - the value of optimal parameters
                    - the associated cost function
        """
        if self.__ready is False:
            if self.verbose:
                print('--> newop Error!')
                print('\tThe problem is not correctly initialized.')
                print('\tRun the "initialize" method before optimizing.')
            return None, None

        # call to SciPy minimizer
        res = minimize(numv, varx0,
                       args=(nproc, self.d_x, self.vdx, self.function),
                       bounds=self.bounds, method='L-BFGS-B',
                       jac=True,
                       options={'maxiter': niter,
                                'maxfun': self.maxfun,
                                'disp': self.verbose}
                       )
        if res.success is True:
            return res.fun, res.x
        else:
            return None, None
