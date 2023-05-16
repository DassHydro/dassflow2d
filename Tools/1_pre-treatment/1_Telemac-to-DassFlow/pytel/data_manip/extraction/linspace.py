"""
Function to lineraize stuff
"""
import numpy as np
import scipy as sc
import scipy.spatial.distance as sc


def linspace_seg(point1, point2, length, end=False):
    """
    Discretize segment into length points

    @param point1 (numpy array) first point of the segment
    @param point2 (numpy array) second point of the segment
    @param length (int) number of discretized points on the segment
    @param end (optional logical) include or not the last point

    @returns (numpy.array)
    """
    v_1 = np.linspace(point1[0], point2[0], length, endpoint=end)
    v_2 = np.linspace(point1[1], point2[1], length, endpoint=end)
    line = np.zeros(shape=[length, 2])
    line[:, 0] = v_1
    line[:, 1] = v_2
    return line


def linspace_poly(poly_points, poly_number):
    """
    Discretize polyline

    @param poly_points (list of numpy array) list of polyline points
    @param poly_number (list of int) number of discretized points for each
    polyline segments

    @returns (numpy.array)
    """
    list_seg = []
    for i in range(len(poly_points)-1):
        lin = linspace_seg(poly_points[i], poly_points[i+1], poly_number[i],
                           end=False)
        list_seg.append(lin)
    length = sum(poly_number)+1
    poly = np.zeros(shape=[length, 2])
    poly[0:poly_number[0], :] = list_seg[0]
    deb = 0
    for i in range(0, len(poly_number)-1):
        deb = deb + poly_number[i]
        poly[deb:(poly_number[i+1]+deb), :] = list_seg[i+1][:]
    poly[-1, :] = poly_points[-1]
    return poly


def curvilinear_abscissa(coord_poly):
    """
    Compute curvilinear abscissa of a polyline

    @param coord_poly (list) coordinates of polyline (list of 2-uple)

    @return curv_absc (list) list of curvilinear abscissa (list of float)
    """
    curv_absc = np.zeros(len(coord_poly))

    for i in range(len(coord_poly)-1):
        coord_point_1 = coord_poly[i, :]
        coord_point_2 = coord_poly[i+1, :]
        curv_absc[i+1] = curv_absc[i] + sc.euclidean(coord_point_1,
                                                     coord_point_2)

    return curv_absc
