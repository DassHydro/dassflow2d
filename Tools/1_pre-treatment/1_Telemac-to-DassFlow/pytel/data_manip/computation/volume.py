#!/usr/bin/env python3
""" Calculating volume difference """

from utils.geometry import get_triangle_area
import numpy as np


def volume_calculation(ikle, variable, coord_x, coord_y, x_ref=0):
    """
    Compute bathymetry volume relatively to zref

    @param ikle (Numpy 2D-array) 2D array of shape (nelem, ndp)
    @param variable (Numpy 1D-array) variable use to compute volume
                                     1D array of size npoin
    @param coord_x (Numpy 1D-array) x coordinate
    @param coord_y (Numpy 1D-array) y coordinate
    @param x_ref (float) volume is calculated relatively to x_ref value
                         (a plan z = x_ref if variable is bottom for example)

    @return the volume, the elemental surface and
                            elemental volume
    """
    surface_point = np.zeros(len(coord_x))
    volume_point = np.zeros(len(coord_x))

    for elem in ikle:
        t_1 = (coord_x[elem[0]], coord_y[elem[0]])
        t_2 = (coord_x[elem[1]], coord_y[elem[1]])
        t_3 = (coord_x[elem[2]], coord_y[elem[2]])

        surface = get_triangle_area(t_1, t_2, t_3)

        surface_point[elem[0]] += surface / 3.0
        surface_point[elem[1]] += surface / 3.0
        surface_point[elem[2]] += surface / 3.0

        volume_point[elem[0]] += surface / 3.0 * (variable[elem[0]] - x_ref)
        volume_point[elem[1]] += surface / 3.0 * (variable[elem[1]] - x_ref)
        volume_point[elem[2]] += surface / 3.0 * (variable[elem[2]] - x_ref)

    return float(volume_point.sum()), surface_point, volume_point


def compute_fv_cell_area(tri):
    """
    Compute finate volume cells area

    @param tri (matplotlib.tri.Triangulation) mesh

    @return area (Numpy 1D-array)
    """
    ikle = tri.triangles
    coord_x = tri.x
    coord_y = tri.y
    area = np.zeros(len(coord_x))

    for elem in ikle:
        t_1 = (coord_x[elem[0]], coord_y[elem[0]])
        t_2 = (coord_x[elem[1]], coord_y[elem[1]])
        t_3 = (coord_x[elem[2]], coord_y[elem[2]])

        surface = get_triangle_area(t_1, t_2, t_3)

        area[elem[0]] += surface / 3.0
        area[elem[1]] += surface / 3.0
        area[elem[2]] += surface / 3.0

    return area
