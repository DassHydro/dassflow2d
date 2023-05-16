#!/usr/bin/env python3
""" Create a triangulation from a grid """

import numpy as np
import matplotlib.tri as mtri


def triangulation_from_data(x, z):
    """
    Extract values of plan in telemac-3d result file for the given variable

    @param x (numpy.array) vector of the curvilinear
    @param z (numpy.array) matrix of elevation (len(x),nplan)

    @returns (matplotlib.tri)

    """
    npoin, nplan = z.shape
    tab_x = np.zeros(nplan*npoin)
    tab_z = np.zeros(nplan*npoin)
    triang = np.zeros((((nplan-1)*(len(x)-1)*2), 3))

    for i, val_x in enumerate(x):
        for j in range(nplan):
            tab_x[i*nplan+j] = val_x
            tab_z[i*nplan+j] = z[i, j]
    for k in range(npoin-1):
        for l in range(nplan-1):
            triang[k*2*(nplan-1)+2*l, 0] = nplan*k+l
            triang[k*2*(nplan-1)+2*l, 1] = nplan*(k+1)+l
            triang[k*2*(nplan-1)+2*l, 2] = nplan*k+(l+1)

            triang[k*2*(nplan-1)+2*l+1, 0] = nplan*k+(l+1)
            triang[k*2*(nplan-1)+2*l+1, 1] = nplan*(k+1)+l
            triang[k*2*(nplan-1)+2*l+1, 2] = nplan*(k+1)+(l+1)

    tri = mtri.Triangulation(tab_x, tab_z, triang)
    return tri
