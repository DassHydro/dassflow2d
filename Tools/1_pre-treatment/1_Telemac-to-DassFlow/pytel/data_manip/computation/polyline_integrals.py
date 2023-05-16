#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

    @brief
    Compute 2D integrals over time across polygonal chains

    @details
    Based on PyTelTools that was originally written by Luc Duron (CNR).
    The extraction has to be done with the dedicated library.
    The integration is performed here.
    Can be used to compute:
    - the area of the wet section along a polygonal chain
    - the flowrate through a polygonal chain
    - the sediment flowrate through a polygonal chain (total load, bedload
      only or suspension load)
    - the flux of any tracer through a polygonal chain
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards other modules
from utils.geometry import get_norm2
from utils.exceptions import TelemacException
import numpy as np

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#


def compute_segments_lengthes(polyline_coords):
    """
    Compute the length of each segment of the polygonal chain,
    and store it in the second point of each segment (the
    first point of the polygonal chain gets a zero length).

    @param polyline_coords: coordinates of the polygonal chain
    """
    length = []
    first_point = None
    for x, y in polyline_coords:
        if first_point is None:  # the first point doesn't have a length
            first_point = (x, y)
            length.append(0)
        else:
            second_point = (x, y)
            length.append(get_norm2(first_point, second_point))
            first_point = second_point
    return length


def compute_segments_normals(polyline_coords):
    """
    Compute the normal to each segment of the polygonal chain
    and store it in the second point of each segment (the
    first point of the polygonal chain gets a zero normal).

    @param polyline_coords: coordinates of the polygonal chain
    """
    normals = []
    prev_x, prev_y = None, None
    for x, y in polyline_coords:
        if prev_x is None:  # the first point doesn't have a normal vector
            prev_x, prev_y = x, y
            normals.append([0, 0])
        else:
            normal_length = get_norm2((x, y), (prev_x, prev_y))
            if normal_length < 10e-10:
                raise TelemacException("The normal "
                                       "length is too small,\
                                        check your mesh and polyline\n")
            normals.append([y-prev_y, prev_x-x]/normal_length)
    return normals


def compute_segments_tangents(polyline_coords):
    """
    Compute the tangents to each segment of the polygonal chain
    and store it in the second point of each segment (the
    first point of the polygonal chain gets a zero normal).

    @param polyline_coords: coordinates of the polygonal chain
    """
    tangents = []
    prev_x, prev_y = None, None
    for x, y in polyline_coords:
        if prev_x is None:  # the first point doesn't have a normal vector
            prev_x, prev_y = x, y
            tangents.append([0, 0])
        else:
            normal_length = get_norm2((prev_y, x), (y, prev_x))
            if normal_length < 10e-10:
                raise TelemacException("The normal "
                                       "length is too small,\
                                        check your mesh and polyline\n")
            tangents.append([x-prev_x, y-prev_y]/normal_length)
    return tangents


def wet_area_2d(polyline_coords, water_depth):
    """
    Compute the wet section over a polygonal chain

    @param polyline_coords: coordinates of the polygonal chain
    @param water_depth: water depth along the polygonal chain
    """
    wet_area = 0
    # first, compute the length of each segment of the polygonal chain
    # and store it in an array
    lengthes = compute_segments_lengthes(polyline_coords)
    # we only loopi from 0 len(lengthes)-1 because the lengthes array has the
    # size of the number of points on the line. Here we loop on the number of
    # segments, which is the number of points - 1
    for i in range(len(lengthes)-1):
        # the extracted values can be nan if the polyline
        # point is not inside the mesh, thus we chose to
        # consider nan values as zero. We make a mistake close to the
        # edge of the domain: for example if water_depth[i] is nan and
        # water_depth[i+1] is not, the segment is crossing the boundary
        # of the domain, but we use the full segment length for the
        # integral calculation while only part of it is actually in the domain
        if np.isnan(water_depth[i][0]):
            water_depth[i] = 0.
        if np.isnan(water_depth[i+1][0]):
            water_depth[i+1] = 0.
        # compute the area of the wet section
        wet_area += (water_depth[i+1] + water_depth[i])*lengthes[i+1]/2.
    return wet_area


def flux_2d(polyline_coords, flux_x, flux_y, scalar=None):
    """
    Compute the 2D flux over a polygonal chain

    @param polyline_coords: coordinates of the polygonal chain
    @param flux_x: value of the flux along x at each point of
        the polygonal chain, it can be
        HU or a solid discharge (QSX, QSBLX or QSSUSPX)
    @param flux_y: value of the flux along y, it can be
        HV or a solid discharge (QSY, QSBLY or QSSUSPY)
    @param scalar: optional, value of a scalar for which
        we want to compute the flux. If it is set this function
        only returns the scalar flux (and not the flow rate or the
        solid discharge)
    """

    flux = 0
    # first, compute the length of each segment of the polygonal chain
    # and store it in an array
    normals = compute_segments_normals(polyline_coords)
    lengthes = compute_segments_lengthes(polyline_coords)
    for i in range(len(lengthes)-1):
        scal_i = 1.
        scal_ip1 = 1.
        # In case flux_2d is called with a scalar argument, fill scal_i
        # and scal_ip1 values. Otherwise they are equal to 1 and do not affect
        # the result
        if scalar is not None:
            scal_i = scalar[i]
            scal_ip1 = scalar[i+1]
        # the extracted values can be nan if the polyline
        # point is not inside the computational mesh, thus we chose to
        # consider nan values as zero. We make a mistake close to the
        # edge of the domain: for example if flux_x[i] is nan and flux_x[i+1]
        # is not, the segment is crossing the boundary of the domain, but we
        # use the full segment length for the integral calculation while only
        # part of it is actually in the domain
        if np.isnan(flux_x[i][0]):
            flux_x[i] = 0.
            flux_y[i] = 0.
            scal_i = 0.
        if np.isnan(flux_x[i+1][0]):
            flux_x[i+1] = 0.
            flux_y[i+1] = 0.
            scal_ip1 = 0.
        # we do not make the check on the scalar and on flux_y, considering
        # that if flux_x is not nan they should be correctly defined

        # compute the fluxes
        product_i = (flux_x[i]*normals[i+1][0]
                     + flux_y[i]*normals[i+1][1]) * scal_i
        product_i_plus_1 = (flux_x[i+1]*normals[i+1][0]
                            + flux_y[i+1]*normals[i+1][1]) * scal_ip1
        flux += (product_i + product_i_plus_1)*lengthes[i+1]/2.
    return flux
