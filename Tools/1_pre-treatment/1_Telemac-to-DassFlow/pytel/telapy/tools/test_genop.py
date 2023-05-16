#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
    Testing the 'genop' package with benchmark problems

    Author(s): Fabrice Zaoui

    Copyright EDF 2016-2018

    Comments : genop... a genetic algorithm optimizer adapted from the original
               algorithm of Scilab module 'optim_ga'
"""

import telapy.tools.genop as genop
import telapy.tools.simul as simul
import unittest
import numpy as np


class TestGenop(unittest.TestCase):
    """ Testing the Genop methods """

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
        mypb = genop.Genop(50)
        mypb.verbose = True
        self.assertEqual(mypb.popsize, 50)
        self.assertEqual(mypb.pbcross, 0.7)
        self.assertEqual(mypb.pbmut, 0.1)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        fcost, _ = mypb.optimize(30, nproc=2)
        print(fcost)
        self.assertTrue(fcost[-1] < 5.e-2)

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
        mypb = genop.Genop(90, 0.69, 0.05)
        mypb.verbose = False
        self.assertEqual(mypb.popsize, 90)
        self.assertEqual(mypb.pbcross, 0.69)
        self.assertEqual(mypb.pbmut, 0.05)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        fcost, _ = mypb.optimize(20)
        self.assertTrue(fcost[-1] < -1.91)

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
        mypb = genop.Genop()
        mypb.verbose = False
        self.assertEqual(mypb.popsize, 25)
        self.assertEqual(mypb.pbcross, 0.7)
        self.assertEqual(mypb.pbmut, 0.1)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        fcost, _ = mypb.optimize(nbgen=100)
        self.assertTrue(np.abs(fcost[-2] - fcost[-1]) < 1.e-12)

    def test_eggholder(self):
        """
        Eggholder function with a large population size
        """
        myfunc = simul.eggholder
        nvar = 2
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -512.
        vbounds[0, 1] = 512.
        vbounds[1, 0] = -512.
        vbounds[1, 1] = 512.
        mypb = genop.Genop(1000)
        mypb.verbose = False
        self.assertEqual(mypb.popsize, 1000)
        self.assertEqual(mypb.pbcross, 0.7)
        self.assertEqual(mypb.pbmut, 0.1)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        fcost, _ = mypb.optimize(20)
        self.assertTrue(fcost[-1] < -900.)

    def test_linear(self):
        """
        Linear function with many iterations
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
        mypb = genop.Genop()
        mypb.verbose = False
        self.assertEqual(mypb.popsize, 25)
        self.assertEqual(mypb.pbcross, 0.7)
        self.assertEqual(mypb.pbmut, 0.1)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        fcost, _ = mypb.optimize(200)
        self.assertTrue(fcost[-1] < -20.)

    def test_abslinear(self):
        """
        Absolute value function with many iterations
        """
        myfunc = simul.abslinear
        nvar = 1
        vbounds = np.zeros((nvar, 2))
        vbounds[0, 0] = -100.
        vbounds[0, 1] = 100.
        mypb = genop.Genop()
        mypb.verbose = False
        self.assertEqual(mypb.popsize, 25)
        self.assertEqual(mypb.pbcross, 0.7)
        self.assertEqual(mypb.pbmut, 0.1)
        error = mypb.initialize(myfunc, nvar, vbounds)
        self.assertEqual(error, 0)
        _, xopt = mypb.optimize(2000)
        self.assertTrue(np.abs(xopt[-1] - 2.) < 1.e-3)


if __name__ == '__main__':
    unittest.main()
