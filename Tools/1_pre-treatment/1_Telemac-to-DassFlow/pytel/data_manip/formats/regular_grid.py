#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium
   @note Functions for regular grid
   @brief
"""
import numpy as np
import matplotlib.tri as mtri
from utils.exceptions import TelemacException


def interpolate_on_grid(tri,
                        data,
                        grid=None,
                        grid_xlim=None,
                        grid_ylim=None,
                        grid_resolution=(500, 500)):
    """
    Interpolation on regular grid (matplotlib.mtri.LinearTriInterpolator)

    @param tri (matplotlib.tri.Triangulation) triangular mesh
    @param data (list(numpy.array) or numpy.array) list of scalar
    data or scalar data
    @param grid (tuple) x and y of interpolation grid
    @param grid_xlim (list(float)) grid limit in x
    @param grid_ylim (list(float)) grid limit in y
    @param grid_resolution (list(float)) number of grid point on x and y
    (default: (500, 500))

    @return interpolated data(list(numpy.array) or numpy.array)
    and grid (tuple)
    """
    if grid is None:
        if grid_xlim is None:
            grid_xlim = [np.min(tri.x), np.max(tri.x)]
        if grid_ylim is None:
            grid_ylim = [np.min(tri.y), np.max(tri.y)]

    if grid is None:
        m_xi, m_yi = np.meshgrid(np.linspace(grid_xlim[0],
                                             grid_xlim[1],
                                             grid_resolution[0]),
                                 np.linspace(grid_ylim[0],
                                             grid_ylim[1],
                                             grid_resolution[1]))
    else:
        m_xi, m_yi = grid

    # Interpolate data x,y on cartesian grid
    if isinstance(data, list):
        data_x_i = []
        for data_x in data:
            data_x_interp = mtri.LinearTriInterpolator(tri, data_x)
            data_x_i.append(data_x_interp(m_xi, m_yi))
    else:
        data_x_interp = mtri.LinearTriInterpolator(tri, data)
        data_x_i = data_x_interp(m_xi, m_yi)

    return data_x_i, (m_xi, m_yi)


def field_diff_on_grid(data1, data2):
    """
    Difference between two fields on the same regulare grid
    @param data1 (list(numpy.array) or numpy.array) list of fields or scalar
    field
    @param data2 (list(numpy.array) or numpy.array) list of fields or field
    @return (list (numpy.array) or numpy.array) difference between two fields
    """
    if isinstance(data1, list) and isinstance(data1, list):
        data_i = []
        for dat1, dat2 in zip(data1, data2):
            data_i.append(np.absolute(dat1-dat2))
    elif not (isinstance(data1, list) and isinstance(dat1, list)):
        data_i = np.absolute(data1 - data2)
    else:
        raise TelemacException('incompatible type between data1 and data2')

    return data_i
