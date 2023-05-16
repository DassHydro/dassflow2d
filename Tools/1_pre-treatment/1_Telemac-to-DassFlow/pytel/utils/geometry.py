r"""
    @note ... this work is based on a collaborative effort
    of the Telemac-Mascaret consortium

    @history 11/11/2011 -- Sebastien E. Bourban
        An accuracy has been introduced because Python does not seem
        to be accurate with sums and multiplications

    @history 07/12/2011 -- Sebastien E. Bourban
        Addition of 3 new geometrical tools:
        + get_segment_line_intersection (different from
          get_segment_intersection)
        + get_plane_equation (of the form Z = a*X + b*Y + c)
        + get_triangle_area

    @history 07/01/2012 -- Sebastien E. Bourban
        Addition of a few geometrical tools, working on angles:
        + get_cone_angle ( based on arctan2 )
        + get_cone_sin_angle ( S = ac.sin(B)/2 = det / 2 )

    @history 14/02/2012 -- Sebastien E. Bourban, Laure C. Grignon
        Addition of the is_inside_poly method to define whether a point is
        inside of a polygon based on the ray casting method

    @brief
        Tools for trivial geometrical operations
"""

# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import math
import numpy as np

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#


def is_ccw(t_1, t_2, t_3):
    """@brief Checks if the element is conterclockwise oriented or not

       @param t_1 (tuple): coordinates of node 1
       @param t_2 (tuple): coordinates of node 2
       @param t_3 (tuple): coordinates of node 3

       @return True if it is oriented counterclockwise, False otherwise
    """
    (x_1, y_1) = t_1
    (x_2, y_2) = t_2
    (x_3, y_3) = t_3

    return (y_3-y_1)*(x_2-x_1) > (y_2-y_1)*(x_3-x_1)


def get_segment_intersection(point_1, point_2, point_3, point_4):
    """
    Returns the coordinate of the point at the intersection
    of two segments, defined by (point_1,point_2)
    and (point_3,point_4) and ratio (exaplined below)

    @param point_1 (2-uple) first point of the first segment
    @param point_2 (2-uple) second point of the first segment
    @param point_3 (2-uple) first point of the second segment
    @param point_4 (2-uple) second point of the second segment

    @return [[x0,y0], ratio] where x0 and y0 are the coordinates of the
    intersect point and ratio that is norm2(point0, point1)/norm2(point1,
    point2)
    """
    (x_1, y_1) = point_1
    (x_2, y_2) = point_2
    (x_3, y_3) = point_3
    (x_4, y_4) = point_4
    det = (x_1-x_2)*(y_3-y_4) - (y_1-y_2)*(x_3-x_4)
    if det == 0:
        return []

    # ~~> Using the mlab tools
   #if mlab.segments_intersect(((x_1, y_1), (x_2, y_2)),
   #                           ((x_3, y_3), (x_4, y_4))):
   #  x_0 = ((x_3-x_4)*(x_1*y_2-y_1*x_2)-(x_1-x_2)*(x_3*y_4-y_3*x_4))/det
   #  y_0 = ((y_3-y_4)*(x_1*y_2-y_1*x_2)-(y_1-y_2)*(x_3*y_4-y_3*x_4))/det
   #  return [[x_0, y_0], get_norm2((x_0, y_0), (x_4, y_4))/\
   #                      get_norm2((x_3, y_3), (x_4, y_4))]
   #return []

    # ~~> Using the clock-wise argument
   #if is_ccw((x_1, y_1), (x_3, y_3), (x_4, y_4)) != \
   #   is_ccw((x_2, y_2), (x_3, y_3), (x_4, y_4)) \
   #  and is_ccw((x_1, y_1), (x_2, y_2), (x_3, y_3)) != \
   #      is_ccw((x_1, y_1), (x_2, y_2), (x_4, y_4)):
   #    x_0 = ((x_3-x_4)*(x_1*y_2-y_1*x_2)-(x_1-x_2)*(x_3*y_4-y_3*x_4))/det
   #    y_0 = ((y_3-y_4)*(x_1*y_2-y_1*x_2)-(y_1-y_2)*(x_3*y_4-y_3*x_4))/det
   #    return [[x_0, y_0], get_norm2((x_0, y_0), (x_4, y_4))/\
   #                        get_norm2((x_3, y_3), (x_4, y_4))]
   #return []

    # ~~> Using the bounding box method
    x_0 = ((x_3-x_4)*(x_1*y_2-y_1*x_2)-(x_1-x_2)*(x_3*y_4-y_3*x_4))/det
    y_0 = ((y_3-y_4)*(x_1*y_2-y_1*x_2)-(y_1-y_2)*(x_3*y_4-y_3*x_4))/det
    accuracy = 0
#   accuracy = np.power(10.0, -5  +np.floor(np.log10(abs(x_1+x_2+x_3+x_4))))
    if (min(x_1, x_2)-x_0) > accuracy or (x_0-max(x_1, x_2)) > accuracy:
        return []
    if (min(x_3, x_4)-x_0) > accuracy or (x_0-max(x_3, x_4)) > accuracy:
        return []
#   accuracy = np.power(10.0, -5  +np.floor(np.log10(abs(y_1+y_2+y_3+y_4))))
    if (min(y_1, y_2)-y_0) > accuracy or (y_0-max(y_1, y_2)) > accuracy:
        return []
    if (min(y_3, y_4)-y_0) > accuracy or (y_0-max(y_3, y_4)) > accuracy:
        return []
    return [[x_0, y_0], get_norm2((x_0, y_0), (x_2, y_2))
            / get_norm2((x_1, y_1), (x_2, y_2))]


def get_segment_line_intersection(t_1, t_2, t_3, t_4):
    """@brief
        Returns the coordinate of the point at the intersection
            of one segments defined by (p_1,p_2) and one line (p_3,p_4)
    """
    (x_1, y_1) = t_1
    (x_2, y_2) = t_2
    (x_3, y_3) = t_3
    (x_4, y_4) = t_4
    det = (x_1-x_2)*(y_3-y_4) - (y_1-y_2)*(x_3-x_4)
    if det == 0:
        return []
    x_0 = ((x_3-x_4)*(x_1*y_2-y_1*x_2)-(x_1-x_2)*(x_3*y_4-y_3*x_4))/det
    y_0 = ((y_3-y_4)*(x_1*y_2-y_1*x_2)-(y_1-y_2)*(x_3*y_4-y_3*x_4))/det
    accuracy = np.power(10.0, -5+np.floor(np.log10(abs(x_1+x_2+x_3+x_4))))
    if (min(x_1, x_2)-x_0) > accuracy or (x_0-max(x_1, x_2)) > accuracy:
        return []
    accuracy = np.power(10.0, -5+np.floor(np.log10(abs(y_1+y_2+y_3+y_4))))
    if (min(y_1, y_2)-y_0) > accuracy or (y_0-max(y_1, y_2)) > accuracy:
        return []

    return [[x_0, y_0]]


def get_norm2(point1, point2):
    """
    Compute norm2 of points given in argument:
    sqrt((y_2-y_1)^2+(x_2-x_1)^2)

    @param point1 (2-uple) Coordinates of the point (x1, y1)
    @param point2 (2-uple) Coordinates of the point (x2, y2)

    @returns (float) The norm
    """
    (x_1, y_1) = point1
    (x_2, y_2) = point2
    return np.sqrt(np.power(x_1-x_2, 2) + np.power(y_1-y_2, 2))


def get_plane_equation(point_1, point_2, point_3):
    """
    Find the equation of the plane defined by 3 points.
    The form of the equation is: Z = a*X + b*Y + c

    @param point_1 (2-uple) coordiantes of the first point
    @param point_2 (2-uple) coordiantes of the second point
    @param point_3 (2-uple) coordiantes of the third point

    @returns (3-uple) (a, b, c)
    """

    (x_1, y_1, z_1) = point_1
    (x_2, y_2, z_2) = point_2
    (x_3, y_3, z_3) = point_3
    det = x_1*(y_2-y_3) + y_1*(x_3-x_2) + (x_2*y_3 - y_2*x_3)
    p_a = (z_1*(y_2-y_3) + z_2*(y_3-y_1) + z_3*(y_1-y_2))/det
    p_b = (z_1*(x_3-x_2) + z_2*(x_1-x_3) + z_3*(x_2-x_1))/det
    p_c = (z_1*(x_2*y_3 - y_2*x_3) + z_2*(y_1*x_3-x_1*y_3) + \
           z_3*(x_1*y_2-y_1*x_2))/det

    return p_a, p_b, p_c


def get_barycentric_weights(point_0, point_1, point_2, point_3):
    """
    Returns the barycentric weights of point point_0 in regard to the three
    points (point_1, point_2, point_3)

    @param point_0 (2-uple) Coordinates of the point for which we want the
    weight
    @param point_1 (2-uple) Coodiantes of the first of the three points
    @param point_2 (2-uple) Coodiantes of the second of the three points
    @param point_3 (2-uple) Coodiantes of the third of the three points

    @returns (3-uple) weight in regard of each point
    """

    (x_1, y_1) = point_1
    (x_2, y_2) = point_2
    (x_3, y_3) = point_3
    (x_o, y_o) = point_0
    det = (y_2-y_3) * (x_1-x_3) - (y_1-y_3) * (x_2-x_3)
    if det == 0.0:
        return 0.0, 0.0, 1.0
    l_1 = ((y_2-y_3) * (x_o-x_3) + (y_o-y_3) * (x_3-x_2))/det
    l_2 = ((y_3-y_1) * (x_o-x_3) + (y_o-y_3) * (x_1-x_3))/det

    return l_1, l_2, 1.0 - l_2 - l_1


def get_distance_point_to_line(point_0, point_1, point_2):
    """
    Returns the distance between a point and a line

    @param point_0 (2-uple) Coodinates of the point
    @param point_1 (2-uple) Coodinates of the first point of the line
    @param point_2 (2-uple) Coodinates of the second point of the line

    @returns (float) the distance

    """

    (x_1, y_1) = point_1
    (x_2, y_2) = point_2
    (x_o, y_o) = point_0

    c_2 = ((x_2-x_1)*(x_2-x_1) + (y_2-y_1)*(y_2-y_1))
    det = (x_2-x_1)*(y_1-y_o) - (x_1-x_o)*(y_2-y_1)

    return abs(det) / math.sqrt(c_2)


def get_triangle_area(point_1, point_2, point_3):
    """
    Compute the area of a triangle

    @param point_1 (2-uple) Coordinates of the first point of the triangle
    @param point_2 (2-uple) Coordinates of the second point of the triangle
    @param point_3 (2-uple) Coordinates of the third point of the triangle

    @returns The area
    """
    # half the vector product
    (x_1, y_1) = point_1
    (x_2, y_2) = point_2
    (x_3, y_3) = point_3
    return 0.5 * abs((x_2-x_1)*(y_3-y_1) - (x_3-x_1)*(y_2-y_1))


def get_cone_sin_angle(t_1, t_2, t_3):
    # TODO: explain what this function does
    (x_1, y_1) = t_1
    (x_2, y_2) = t_2
    (x_3, y_3) = t_3
    # S = ac.sin(B)/2 = det / 2
    a_2 = ((x_2-x_3)*(x_2-x_3) + (y_2-y_3)*(y_2-y_3))
    c_2 = ((x_2-x_1)*(x_2-x_1) + (y_2-y_1)*(y_2-y_1))
    a_c = np.sqrt(a_2*c_2)
    det = (x_1-x_2)*(y_3-y_2) - (x_3-x_2)*(y_1-y_2)  # A to C
    return det / a_c


def get_cone_angle(t_1, t_2, t_3):
    # TODO: explain what this function does
    (x_1, y_1) = t_1
    (x_2, y_2) = t_2
    (x_3, y_3) = t_3
    return np.arctan2(y_2-y_3, x_2-x_3) - np.arctan2(y_2-y_1, x_2-x_1)


def is_inside_triangle(t_o, t_1, t_2, t_3,
                       size=5, nomatter=False):
    # TODO: This function returns an array but is suppose to say if a point is
    # inside a triangle ????
    (x_1, y_1) = t_1
    (x_2, y_2) = t_2
    (x_3, y_3) = t_3
    (x_o, y_o) = t_o
    # ~~> Taking sides
    l_1 = (x_o-x_2)*(y_1-y_2)-(x_1-x_2)*(y_o-y_2) < 0.0
    l_2 = (x_o-x_3)*(y_2-y_3)-(x_2-x_3)*(y_o-y_3) < 0.0
    l_3 = (x_o-x_1)*(y_3-y_1)-(x_3-x_1)*(y_o-y_1) < 0.0
    if l_1 == l_2 == l_3:
        return get_barycentric_weights((x_o, y_o), (x_1, y_1), (x_2, y_2),
                                       (x_3, y_3))
    if nomatter:
        return get_barycentric_weights((x_o, y_o), (x_1, y_1), (x_2, y_2),
                                       (x_3, y_3))
    return []

    # ~~> Using barycentric weight
   #l_1 ,l_2, l_3 = get_barycentric_weights((x_o, y_o), (x_1, y_1),
   #                                        (x_2, y_2), (x_3, y_3))
   #accuracy = np.power(10.0, -size+np.floor(np.log10(abs(l_1+l_2+l_3))))
   #if l_1 >= -accuracy and l_1 <= 1.0+accuracy and \
   #  l_2 >= -accuracy and l_2 <= 1.0+accuracy and \
   #  l_3 >= -accuracy and l_3 <= 1.0+accuracy:
   #     return [l_1, l_2, l_3]
   #if nomatter:
   #  return [l_1, l_2, l_3]
   #return []


def is_inside_poly(t_o, poly, close=True):
    """
    Check if a point is in a polygon using Ray Casting Method

    @param t_o (2-uple) Coordinates of the point
    @param poly (list) list of point of the polygon
    @param close (boolean) if true also check if the point is near the boundary
    of the polygon

    @returns (boolean) True if inside
    """
    # by the "Ray Casting Method"
    (x_o, y_o) = t_o
    inside = False
    p_1 = poly[0]
    for j in range(len(poly)+1):
        p_2 = poly[j % len(poly)]
        if y_o >= min(p_1[1], p_2[1]):
            if y_o <= max(p_1[1], p_2[1]):
                if x_o <= max(p_1[0], p_2[0]):
                    if p_1[1] != p_2[1]:
                        xints = (y_o-p_1[1])*(p_2[0]-p_1[0])/(p_2[1]-p_1[1])+\
                                p_1[0]
                    if p_1[0] == p_2[0] or x_o <= xints:
                        inside = not inside
        p_1 = p_2
    if close:
        for p_1 in poly:
            if is_close([x_o, y_o], p_1, size=10):
                inside = True
    return inside


def is_close(p_1, p_2, size=5):
    """
    Checks if two points are near each other

    @param p_1 (2-uple) Coordianate of the first point
    @param p_2 (2-uple) Coordianate of the second point
    @param size (int) Additional offset for accuracy (10^-size)
    """

    if p_2 == [] or p_1 == []:
        return False
    tmp = 1.e-5 + abs(max(p_1)+max(p_2))
    accuracy = np.power(10.0, -size+np.floor(np.log10(tmp)))

    return get_norm2(p_1[0:2], p_2[0:2]) < accuracy

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#


__author__ = "Sebastien E. Bourban"
__date__ = "$15-Nov-2011 08:51:29$"
