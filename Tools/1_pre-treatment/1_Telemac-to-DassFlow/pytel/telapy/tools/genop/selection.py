# -*- coding: utf-8 -*-
"""
Selection operator

Author(s) : Fabrice Zaoui

Copyright EDF 2016-2017

:param 'pop': the genetic population
:param 'fobjpop': the cost functions
:param 'nbcouples' : the number of couples to set
:param 'efficiency': selection parameter
:return: two arrays of coupled individuals with their cost function
"""

import numpy as np


def selection(pop, fobjpop, nbcouples, efficiency):
    """
    Selection genetic operator
    """
    (_, nvar) = pop.shape
    wheel = np.cumsum(efficiency)
    func1 = np.array([]).reshape(0, 1)
    func2 = np.array([]).reshape(0, 1)
    ind1 = np.array([]).reshape(0, nvar)
    ind2 = np.array([]).reshape(0, nvar)
    for _ in range(0, nbcouples):
        # selection of the first individual in the couple
        shoot = np.random.uniform(0, wheel[-1], (1, 1))
        index = np.nonzero(shoot <= wheel)
        ind1 = np.vstack([ind1, pop[index[1][0]]])
        func1 = np.vstack([func1, fobjpop[index[1][0]]])
        # selection of the second individual in the couple
        shoot = np.random.uniform(0, wheel[-1], (1, 1))
        index = np.nonzero(shoot <= wheel)
        ind2 = np.vstack([ind2, pop[index[1][0]]])
        func2 = np.vstack([func2, fobjpop[index[1][0]]])
    return ind1, func1, ind2, func2


def elitist(pop, ind1, ind2, fvalpop, func1, func2, pressure):
    """
    Elitist improvment
    """
    (popsize, _) = pop.shape
    pop = np.vstack((pop, ind1))
    pop = np.vstack((pop, ind2))
    totalpop = pop
    fvalpop = np.vstack((fvalpop, func1))
    fvalpop = np.vstack((fvalpop, func2))
    totalfobj = fvalpop
    # normalization of the efficiency
    fmin = np.min(totalfobj)
    fmax = np.max(totalfobj)
    efficiency = (1. - pressure) * (fmax - totalfobj) / \
        np.ndarray.max(np.array([fmax - fmin, np.spacing(1)]), axis=None) + \
        pressure
    order = np.argsort(efficiency, axis=0)
    order = np.flipud(order)
    order = order.reshape((order.size,))
    efficiency = np.array(efficiency)[order]
    efficiency = efficiency[0:popsize]
    # extraction and selection of the phenotype
    totalfobj = np.array(totalfobj)[order]
    fvalpop = totalfobj[0:popsize]
    # extraction and selection of the genotype
    totalpop = np.array(totalpop)[order]
    pop = totalpop[0:popsize]
    return pop, fvalpop, efficiency
