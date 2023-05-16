# -*- coding: utf-8 -*-
"""
Evaluate the cost function and gradients

Author(s) : Fabrice Zaoui

Copyright EDF 2017-2018

"""

import multiprocessing as mp
import numpy as np


def numv(x, nproc, d_x, vdx, fname):
    """
    Compute the cost function and 1st order derivative values

    :param 'x': the local point
    :param 'nproc': the number of processors to use
    :param 'd_x', 'vdx': FD step sizes
    :param 'fname': name of the cost function

    :return: the cost function and its first derivatives
    """
    dimx = x.size
    x_p = np.zeros((dimx, 1))
    if vdx is None:
        for i in range(dimx):
            x_p[i] = x[i] + d_x
    else:
        for i in range(dimx):
            x_p[i] = x[i] + vdx[i]
    pop = np.array([]).reshape(0, dimx)
    pop = np.vstack([pop, x])
    for i in range(dimx):
        pop = np.vstack([pop, x])
        pop[1+i, i] = x_p[i]
    pool = mp.Pool(processes=nproc)
    feval = pool.map(fname, pop)
    pool.close()
    pool.join()
    val = np.asarray(feval[0])
    if vdx is None:
        jac = (np.asarray(feval[1:]) - val) / d_x
    else:
        jac = np.zeros(dimx)
        for i in range(dimx):
            jac[i] = (np.asarray(feval[i+1]) - val) / vdx[i]
    return val, jac
