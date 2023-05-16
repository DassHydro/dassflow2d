#  Copyright (C) 2012-2015 EDF
#
#  This file is part of SALOME HYDRO module.
#
#  SALOME HYDRO module is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SALOME HYDRO module is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with SALOME HYDRO module.  If not, see <http://www.gnu.org/licenses/>.
"""
   Moduel handling polygon
"""
import math

EPS = 10**-9


def is_in_polygon(pt_x, pt_y, poly):
    """
    Returns True if the point (x,y) is in the polygon poly
    Warning the polygon must be distinct (no indentical points)

    param pt_x X coordinate of the point
    param pt_y Y coordinate of the point
    param poly List of the points representing the polygon
    """
    angle = -1.0
    while True:
        angle += 1.0
        if angle > 360.0:
            # Special case of a point on the contour
            return True
        dist_a = math.cos(angle*math.pi/180.0)
        dist_b = math.sin(angle*math.pi/180.0)
        nsect = 0
        # Loop on all the segment of the polygon
        poly_len = len(poly)
        next_angle = False
        for i in range(0, poly_len):
            xdep, ydep = poly[i]
            # next element
            xarr, yarr = poly[(i+1) % poly_len]
            # Case the point is on the polygon
            if (abs(pt_x-xdep) < EPS) and (abs(pt_y-ydep) < EPS):
                nsect = 1
                break
            det = dist_a*(ydep-yarr) - dist_b*(xdep-xarr)
            if abs(det) < EPS:
                next_angle = True
                break

            mu_val = ((xdep-pt_x)*(ydep-yarr) - (ydep-pt_y)*(xdep-xarr)) / det
            llambda = (dist_a*(ydep-pt_y) - dist_b*(xdep-pt_x)) / det
            # if the intersection point is a vertex, increases the angle
            # otherwise the point would be counted twice instead of just once
            if(abs(pt_x+dist_a*mu_val-xdep) <= EPS
               and abs(pt_y+dist_b*mu_val-ydep) <= EPS) or \
              (abs(pt_x+dist_a*mu_val-xarr) <= EPS
               and abs(pt_y+dist_b*mu_val-yarr) <= EPS):
                next_angle = True
                break
            if mu_val >= -EPS and -EPS <= llambda <= (1.0+EPS):
                nsect += 1
        if not next_angle:
            break

    return (nsect % 2) == 1


def import_poly_from_file(filename, sep=','):
    """
    import the polygon coordinates from file

    param filename File which contains the polygon coordinates
    param poly List of the points representing the polygon
    """
    # Building polygon from file
    poly = []
    with open(filename, 'r') as fle:
        for line in fle:
            x_val, y_val = line.split(sep)
            poly.append((float(x_val), float(y_val)))
    return poly

def points_in_poly(points, polyline):
    """
    Test if points are in polyline using matplotlib path

    @param points (np.array of shape (n,2)) contains x,y coordonates of points
    @param polyline (np.array of shape (m,2)) contains x,y coordonates of
    points in the polyline
    @ return inside_poly (list of bool of shape (n)) true if point is in
    polyline
    """
    from matplotlib.path import Path
    poly_path = Path(polyline)
    inside_poly = poly_path.contains_points(points)
    return inside_poly


if __name__ == "__main__":
    MY_POLY = [(0.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 0.0)]
    print(repr(is_in_polygon(1.0, 1.0, MY_POLY)) + ' True')
    print(repr(is_in_polygon(0.0, 0.0, MY_POLY)) + ' True')
    print(repr(is_in_polygon(-1.0, 1.0, MY_POLY)) + ' False')
    print(repr(is_in_polygon(2.0, 0.0, MY_POLY)) + ' True')
