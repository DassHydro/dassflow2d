#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
    Testing the 'newop' package with benchmark problems

    Author(s): Fabrice Zaoui

    Copyright EDF 2017-2018

    Comments : A quasi-Newton optimizer based on the SciPy 'minimizer' function
"""

import unittest
import telapy.tools.newop as newop
import telapy.tools.simul as simul
import numpy as np


class TestNewop(unittest.TestCase):
    """ Testing the Newop methods """

    def test_rosenbrock(self):
        """
        Rosenbrock function with verbose mode and parallelism
        """
        myfunc = simul.rosenbrock
        nvar = 2
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -5.
        vbounds[0, 1] = 5.
        vbounds[1, 0] = -5.
        vbounds[1, 1] = 5.
        mypb = newop.Newop()
        mypb.verbose = True
        self.assertEqual(mypb.d_x, 1.e-6)
        self.assertEqual(mypb.maxfun, 2000)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        varx0 = (vbounds[:, 0] + vbounds[:, 1]) / 2.
        fcost, _ = mypb.optimize(varx0, 50, nproc=2)
        self.assertTrue(fcost < 1.e-6)

    def test_maccormick(self):
        """
        Mc Cormick function with non default algorithm parameters
        """
        myfunc = simul.maccormick
        nvar = 2
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -1.5
        vbounds[0, 1] = 4.
        vbounds[1, 0] = -3.
        vbounds[1, 1] = 4.
        mypb = newop.Newop()
        mypb.verbose = False
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        varx0 = (vbounds[:, 0] + vbounds[:, 1]) / 2.
        fcost, _ = mypb.optimize(varx0, 20)
        self.assertTrue(fcost < -1.91)

    def test_easom(self):
        """
        Easom function
        """
        myfunc = simul.easom
        nvar = 2
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -100.
        vbounds[0, 1] = 100.
        vbounds[1, 0] = -100.
        vbounds[1, 1] = 100.
        mypb = newop.Newop()
        mypb.verbose = False
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        varx0 = (vbounds[:, 0] + vbounds[:, 1]) / 2.
        varx0[0] = 2.5
        varx0[1] = 2.5
        fcost, _ = mypb.optimize(varx0, niter=5)
        self.assertTrue(np.abs(fcost + 1.) < 1.e-12)

    def test_eggholder(self):
        """
        Eggholder function with two steps for the numerical derivatives
        """
        myfunc = simul.eggholder
        nvar = 2
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -512.
        vbounds[0, 1] = 512.
        vbounds[1, 0] = -512.
        vbounds[1, 1] = 512.
        mypb = newop.Newop()
        mypb.verbose = False
        vdx = np.array([1.e-6, 1.e-8])
        error = mypb.initialize(myfunc, nvar, vbounds, vdx)
        self.assertEqual(error, 0)
        varx0 = (vbounds[:, 0] + vbounds[:, 1]) / 2.
        fcost, _ = mypb.optimize(varx0, 10)
        self.assertTrue(fcost < -66.)

    def test_linear(self):
        """
        Linear function
        """
        myfunc = simul.linear
        nvar = 3
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -3.14
        vbounds[0, 1] = 3.14
        vbounds[1, 0] = 0.
        vbounds[1, 1] = 5.
        vbounds[2, 0] = -10.
        vbounds[2, 1] = 20.
        mypb = newop.Newop()
        mypb.verbose = False
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        varx0 = (vbounds[:, 0] + vbounds[:, 1]) / 2.
        fcost, _ = mypb.optimize(varx0, 10)
        self.assertTrue(fcost < -23.1)

    def test_abslinear(self):
        """
        Absolute value function (non-differentiable)
        """
        myfunc = simul.abslinear
        nvar = 1
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -100.
        vbounds[0, 1] = 100.
        mypb = newop.Newop()
        mypb.verbose = False
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        varx0 = (vbounds[:, 0] + vbounds[:, 1]) / 2.
        fcost, _ = mypb.optimize(varx0)
        self.assertTrue(fcost is None)


if __name__ == '__main__':
    unittest.main()
