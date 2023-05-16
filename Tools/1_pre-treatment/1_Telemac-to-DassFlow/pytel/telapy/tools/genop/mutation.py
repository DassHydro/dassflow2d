# -*- coding: utf-8 -*-
"""
Mutation operator

Author(s) : Fabrice Zaoui

Copyright EDF 2016-2017

:param 'ind1': first individual of couples
:param 'ind2': second individual of couples
:param 'todo1': consider the re-computation of the cost function for 'ind1'
:param 'todo2': consider the re-computation of the cost function for 'ind2'
:param 'bounds' : lower and upper bounds of variables
:param 'pbmut': mutation probability
:return: updated populations ind1 and ind2
"""

import numpy as np


def mutation(ind1, ind2, todo1, todo2, bounds, pbmut):
    """
    Mutation genetic operator
    """
    delta = 0.1
    (nbcouples, nvar) = ind1.shape
    for i in range(0, nbcouples):
        if pbmut > np.random.uniform(0, 1, (1, 1)):
            xvar1 = ind1[i, :]
            mut1 = xvar1 + 2. * delta * np.random.uniform(0, 1, (1, nvar)) - \
                delta * np.ones((1, nvar))
            mut1 = np.maximum(np.minimum(mut1, bounds[:, 1]), bounds[:, 0])
            ind1[i, :] = mut1
            todo1[i] = 1
        if pbmut > np.random.uniform(0, 1, (1, 1)):
            xvar2 = ind2[i, :]
            mut2 = xvar2 + 2. * delta * np.random.uniform(0, 1, (1, nvar)) - \
                delta * np.ones((1, nvar))
            mut2 = np.maximum(np.minimum(mut2, bounds[:, 1]), bounds[:, 0])
            ind2[i, :] = mut2
            todo2[i] = 1

    return ind1, ind2, todo1, todo2
