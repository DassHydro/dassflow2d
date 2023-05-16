# -*- coding: utf-8 -*-
"""
Crossover operator

Author(s) : Fabrice Zaoui

Copyright EDF 2016-2017

:param 'ind1': first individual of couples
:param 'ind2': second individual of couples
:param 'bounds' : lower and upper bounds of variables
:param 'pbcross': crossover probability
:return: updated populations ind1 and ind2
"""

import numpy as np


def crossover(ind1, ind2, bounds, pbcross):
    """
    Crossover genetic operator
    """
    (nbcouples, nvar) = ind1.shape
    todo1 = np.zeros((nbcouples, 1))
    todo2 = np.zeros((nbcouples, 1))
    for i in range(0, nbcouples):
        if pbcross > np.random.uniform(0, 1, (1, 1)):
            mix = np.random.uniform(0, 1, (1, 1))
            cross1 = mix * ind1[i, :] + (1. - mix) * ind2[i, :]
            cross2 = (1. - mix) * ind1[i, :] + mix * ind2[i, :]
            cross1 = np.maximum(np.minimum(cross1, bounds[:, 1]), bounds[:, 0])
            cross2 = np.maximum(np.minimum(cross2, bounds[:, 1]), bounds[:, 0])
            ind1[i, :] = cross1
            ind2[i, :] = cross2
            todo1[i] = 1
            todo2[i] = 1

    return ind1, ind2, todo1, todo2
