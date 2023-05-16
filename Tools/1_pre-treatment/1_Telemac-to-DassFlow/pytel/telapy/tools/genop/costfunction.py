# -*- coding: utf-8 -*-
"""
Evaluate the cost functions

Author(s) : Fabrice Zaoui

Copyright EDF 2016-2017

:param 'pop': the genetic population
:param 'npop': the size of the population
:param 'nvar': the number of variables
:param 'fname': name of the cost function
:param 'ind1': first individual of couples
:param 'ind2': second individual of couples
:param 'func1': cost function for the first individual
:param 'func2': cost function for the second individual
:param 'todo1': consider the re-computation of the cost function for 'ind1'
:param 'todo2': consider the re-computation of the cost function for 'ind2'
:return: cost function values
"""

import multiprocessing as mp
import numpy as np


def cost(pop, npop, nvar, fname, nproc):
    """
    Compute the cost function value for each member of pop

    :param 'pop': the genetic population
    :param 'npop': the size of the population
    :param 'nvar': the number of variables
    :param 'fname': name of the cost function
    :param 'nproc': Number of process

    :return: cost function values
    """
    # serial computations
    if nproc == 1:
        ncalls = 0
        feval = np.array([]).reshape(0, 1)
        for i in range(0, npop):
            res = fname(pop[i])
            feval = np.vstack([feval, res])
            ncalls = ncalls + 1
        return feval, ncalls
    # parallel computations
    else:
        pool = mp.Pool(processes=nproc)
        feval = pool.map(fname, pop)
        pool.close()
        pool.join()
        ncalls = npop
        return np.asarray(feval).reshape(npop, 1), ncalls


def updatecost(ind1, ind2, func1, func2, todo1, todo2, fname, nproc):
    """
    Updating some cost function values pointed by the 'todo' lists

    :param 'ind1': first individual of couples
    :param 'ind2': second individual of couples
    :param 'func1': cost function for the first individual
    :param 'func2': cost function for the second individual
    :param 'todo1': consider the re-computation of the cost function for 'ind1'
    :param 'todo2': consider the re-computation of the cost function for 'ind2'
    :param 'fname': name of the cost function
    :param 'nproc': Number of process

    :return: cost function values
    """
    # serial computations
    if nproc == 1:
        ncalls = 0
        (nbcouples, nvar) = ind1.shape
        for i in range(0, nbcouples):
            if todo1[i]:
                func1[i] = fname(ind1[i])
                ncalls = ncalls + 1
            if todo2[i]:
                func2[i] = fname(ind2[i])
                ncalls = ncalls + 1
    # parallel computations
    else:
        ncalls = 0
        (nbcouples, nvar) = ind1.shape
        ind1tocompute = np.array([]).reshape(0, nvar)
        ind2tocompute = np.array([]).reshape(0, nvar)
        for i in range(0, nbcouples):
            if todo1[i]:
                ind1tocompute = np.vstack([ind1tocompute, ind1[i]])
                ncalls = ncalls + 1
            if todo2[i]:
                ind2tocompute = np.vstack([ind2tocompute, ind2[i]])
                ncalls = ncalls + 1
        pool = mp.Pool(processes=nproc)
        gfunc1 = pool.map(fname, ind1tocompute)
        gfunc2 = pool.map(fname, ind2tocompute)
        pool.close()
        pool.join()
        gfunc1 = np.asarray(gfunc1)
        gfunc2 = np.asarray(gfunc2)
        kcount1 = 0
        kcount2 = 0
        for i in range(0, nbcouples):
            if todo1[i]:
                func1[i] = gfunc1[kcount1]
                kcount1 = kcount1 + 1
            if todo2[i]:
                func2[i] = gfunc2[kcount2]
                kcount2 = kcount2 + 1
    return func1, func2, ncalls
