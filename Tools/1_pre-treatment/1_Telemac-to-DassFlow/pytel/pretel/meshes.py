r"""@author Sebastien E. Bourban

    @note ... this work is based on a collaborative effort between
  .________.                                                          ,--.
  |        |                                                      .  (  (
  |,-.    /   HR Wallingford                EDF - LNHE           / \_ \_/ .--.
  /   \  /    Howbery Park,                 6, quai Watier       \   )   /_   )
    ,.  `'     Wallingford, Oxfordshire      78401 Cedex           `-'_  __ `--
  /  \   /    OX10 8BA, United Kingdom      Chatou, France        __/ \ \ `.
 /    `-'|    www.hrwallingford.com         innovation.edf.com   |    )  )  )
!________!                                                        `--'   `--

    @history 12/12/2012 -- Sebastien E. Bourban
      Many methods developped for application to meshes. The latest
          one being about subdivision of meshes.

    @brief
        Tools for sampling and interpolating through triangular meshes

    @details
            Contains ...

    @history 20/06/2013 -- Sebastien E. Bourban
        A new method, slice_mesh, now replaces cross_mesh and all of the
            Ray Tracing algorithms. The later will remain for the fame and
            maybe for future uses, but slive_mesh should now be used.

    @history 10/11/2015 -- Juliette Parisi and Sebastien E. Bourban
        A new method, tessellate_poly, has now been implmented.
        It returns a MESH (including the ikle and ipobo) of the waters within
            the largest polygon, having taken out all islands.
        TODO:
         - It should be capable of taking in open lines to constrain the
         mesh to these.
         - It should be parallelised fairly easily.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
import math
from scipy.spatial import cKDTree
import matplotlib.path as mplPath
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation
# ~~> dependencies towards other modules
# ~~> dependencies towards other pytel/modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.geometry import is_ccw, get_segment_intersection, \
                           get_barycentric_weights, is_inside_triangle, \
                           get_distance_point_to_line
from utils.exceptions import TelemacException
from pretel.polygons import is_clockwise, join_segments

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

def near_locate_mesh(xyo, ikle, meshx, meshy, tree=None):
    """
    Requires the scipy.spatial and the matplotlib.tri packages to be loaded.
     - Will use already computed tree or re-create it if necessary.
     - Will use already computed neighbourhood or re-create it if necessary.
    This function return the element number for the triangle including
     xyo=(xo,yo) or -1 if the (xo,yo) is outside the mesh
    Return: the element, the barycentric weights, and the tree and the
              neighbourhood if computed
    """
    # ~~> Create the KDTree of the iso-barycentres
    if tree == None:
        isoxy = np.column_stack((np.sum(meshx[ikle], axis=1)/3.0,
                                 np.sum(meshy[ikle], axis=1)/3.0))
        tree = cKDTree(isoxy)
    # ~~> Find the indices corresponding to the nearest elements to the points
    inear = -1
    for my_d, i in zip(*tree.query(xyo, 8)):
        a_x, b_x, c_x = meshx[ikle[i]]
        a_y, b_y, c_y = meshy[ikle[i]]
        w = is_inside_triangle(xyo, (a_x, a_y), (b_x, b_y), (c_x, c_y),
                               nomatter=True)
        if w != []:
            return i, w, tree
        if inear < 0:
            inear = i
            dnear = my_d
        if dnear > my_d:
            inear = i
            dnear = my_d

    # ~~> Find the indices and weights corresponding to the element containing
    # the point
    a_x, b_x, c_x = meshx[ikle[inear]]
    a_y, b_y, c_y = meshy[ikle[inear]]

    return inear, is_inside_triangle(xyo, (a_x, a_y), (b_x, b_y), (c_x, c_y),
                                     nomatter=False), tree

def dicho_locate_mesh(rank, e_1, xy1, e_2, xy2, ikle, meshx, meshy, tree):
    """
    Will find at least one point between xy1 and xy2 that is within the mesh
    """
    # ~~ Position the middle point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xyo = [(xy1[0]+xy2[0])/2.0, (xy1[1]+xy2[1])/2.0]
    e_o, b_o, tree = near_locate_mesh(xyo, ikle, meshx, meshy, tree)
    if b_o != []:
        return True, e_o, xyo, b_o

    # ~~ Limit the number of useless dichotomies ~~~~~~~~~~~~~~~~~~~~
    rank = rank + 1
    if rank > 3:
        return False, e_o, xyo, b_o

    # ~~ Sub-segments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    found, e_j, xyj, b_j = dicho_locate_mesh(rank, e_1, xy1, e_o, xyo,
                                             ikle, meshx, meshy, tree)
    if found:
        return found, e_j, xyj, b_j
    found, e_j, xyj, b_j = dicho_locate_mesh(rank, e_o, xyo, e_2, xy2,
                                             ikle, meshx, meshy, tree)
    if found:
        return found, e_j, xyj, b_j

    return False, e_o, xyo, b_o

def xy_trace_mesh(inear, xyi, xyo, ikle, meshx, meshy, neighbours=None):
    """
    Requires the matplotlib.tri package to be loaded.
     - Will use already computed neighbourhood or re-create it if necessary.
    This function return the element number for the triangle including
    xyo=(xo,yo) or -1 if the (xo,yo) is outside the mesh. It creates the
    neighbouring connectivity map and ray-traces from xyi to xyo
    Return: whether xyo was found within the mesh, the (nearest) element, the
        associated barycentric weights, and the neighbourhood if computed
    """
    if neighbours == None:
        neighbours = Triangulation(meshx, meshy, ikle)\
                     .get_cpp_triangulation().get_neighbors()
    found, ray = trace_ray2xy(ikle, meshx, meshy, neighbours, inear, xyi,
                              inear, xyo)

    return found, ray, neighbours

def subdivide_mesh4(ikle, meshx, meshy):
    """
    Requires the matplotlib.tri package to be loaded.
     - Will use already computed edges or re-create it if necessary.
    This function return a new tuple ikle,meshx,meshy where each triangle has
    been subdivided in 4.
    Note: in order to search for segments (split node number on segment), these
        are ordered b_y the minimum of their node number. Thus two neighbouring
        elements see the common segment the same way, from the minimum number to
        the maximum number.
    """
    # ~~> Singling out edges
    edges = Triangulation(meshx, meshy, ikle).get_cpp_triangulation()\
                                             .get_edges()

    # ~~> Memory allocation for new MESH
    ielem = len(ikle)
    ipoin = len(meshx)
    iedge = len(edges)
    # you subdivide every elements b_y 4
    jkle = np.zeros((ielem*4, 3), dtype=np.int)
    # you add one point on every edges
    meshj = np.zeros((iedge, 2), dtype=np.int)
    # n_o 3-node interpolation not necessary
    meshk = np.array([], dtype=np.int)

    # ~~> Lookup tables for node numbering on common edges
    p_a, p_b = edges.T
    k1b, k1a = np.sort(np.take(ikle, [0, 1], axis=1)).T
    indx1 = np.searchsorted(p_a, k1a)
    jndx1 = np.searchsorted(p_a, k1a, side='right')
    k2b, k2a = np.sort(np.take(ikle, [1, 2], axis=1)).T
    indx2 = np.searchsorted(p_a, k2a)
    jndx2 = np.searchsorted(p_a, k2a, side='right')
    k3b, k3a = np.sort(np.take(ikle, [2, 0], axis=1)).T
    indx3 = np.searchsorted(p_a, k3a)
    jndx3 = np.searchsorted(p_a, k3a, side='right')

    # ~~> Building one triangle at a time /!\ Please get this loop parallelised
    j = 0
    for i in range(ielem):
        k_1 = indx1[i]+np.searchsorted(p_b[indx1[i]:jndx1[i]], k1b[i])
        k_2 = indx2[i]+np.searchsorted(p_b[indx2[i]:jndx2[i]], k2b[i])
        k_3 = indx3[i]+np.searchsorted(p_b[indx3[i]:jndx3[i]], k3b[i])
        # ~~> New connectivity jkle
        jkle[j] = [ikle[i][0], ipoin+k_1, ipoin+k_3]
        jkle[j+1] = [ikle[i][1], ipoin+k_2, ipoin+k_1]
        jkle[j+2] = [ikle[i][2], ipoin+k_3, ipoin+k_2]
        jkle[j+3] = [ipoin+k_1, ipoin+k_2, ipoin+k_3]
        # ~~> New interpolation references for values and coordinates
        meshj[k_1] = [ikle[i][0], ikle[i][1]]
        meshj[k_2] = [ikle[i][1], ikle[i][2]]
        meshj[k_3] = [ikle[i][2], ikle[i][0]]
        j += 4

    # ~~> Reset ipobo while you are at it
    meshx = np.resize(meshx, ipoin+iedge)
    meshy = np.resize(meshy, ipoin+iedge)
    meshx[ipoin:] = np.sum(meshx[meshj], axis=1)/2.
    meshy[ipoin:] = np.sum(meshy[meshj], axis=1)/2.
    neighbours = Triangulation(meshx, meshy, jkle).get_cpp_triangulation()\
                                                  .get_neighbors()
    jpobo = np.zeros(ipoin+iedge, np.int)
    for n in range(ielem*4):
        s_1, s_2, s_3 = neighbours[n]
        e_1, e_2, e_3 = jkle[n]
        if s_1 < 0:
            jpobo[e_1] = e_1+1
            jpobo[e_2] = e_2+1
        if s_2 < 0:
            jpobo[e_2] = e_2+1
            jpobo[e_3] = e_3+1
        if s_3 < 0:
            jpobo[e_3] = e_3+1
            jpobo[e_1] = e_1+1

    return jkle, meshx, meshy, jpobo, meshj, meshk

def subdivide_mesh3(ikle, meshx, meshy):
    """
    Requires the matplotlib.tri package to be loaded.
     - Will use already computed edges or re-create it if necessary.
    This function return a new tuple ikle,meshx,meshy where each triangle has
    been subdivided in 3, and where possible re-combined with neighbouring
    triangles.
    """
    # ~~> Singling out edges
    triangles = Triangulation(meshx, meshy, ikle).get_cpp_triangulation()
    edges = triangles.get_edges()
    neighbours = triangles.get_neighbors()

    # ~~> Memory allocation for new MESH
    ielem = len(ikle)
    ipoin = len(meshx)
    _ = len(edges)
    isbnd = 3*ielem - np.count_nonzero(neighbours+1)
    # you subdivide every elements b_y 3 swap + boundary edges
    jkle = np.zeros((ielem*3+isbnd, 3), dtype=np.int)
    # you add one point in the middle of every triangle
    meshk = np.zeros((ielem, 3), dtype=np.int)
    # you add one point on every boundary edges
    meshj = np.zeros((isbnd, 2), dtype=np.int)

    # ~~> Lookup tables for node numbering on common edges
    p_a, _ = edges.T
    _, k1a = np.sort(np.take(ikle, [0, 1], axis=1)).T
    _ = np.searchsorted(p_a, k1a)
    _ = np.searchsorted(p_a, k1a, side='right')
    _, k2a = np.sort(np.take(ikle, [1, 2], axis=1)).T
    _ = np.searchsorted(p_a, k2a)
    _ = np.searchsorted(p_a, k2a, side='right')
    _, k3a = np.sort(np.take(ikle, [2, 0], axis=1)).T
    _ = np.searchsorted(p_a, k3a)
    _ = np.searchsorted(p_a, k3a, side='right')

    # ~~> Building one triangle at a time /!\ Please get this loop parallelised
    j = 0
    k = 0
    for i in range(ielem):
        s_1, s_2, s_3 = neighbours[i]
        meshk[i] = ikle[i]
        # ~~> New boundary segments
        if s_1 < 0:
            jkle[j] = [ikle[i][0], ipoin+ielem+k, ipoin+i]
            jkle[3*ielem+k] = [ipoin+i, ipoin+ielem+k, ikle[i][1]]
            meshj[k] = [ikle[i][0], ikle[i][1]]
            j += 1
            k += 1
        # ~~> New inside segment, swapped
        else:
            # look for the connection with s_1
            e_1, e_2, e_3 = neighbours[s_1]
            n = -1
            if e_1 == i:
                n = 0
            if e_2 == i:
                n = 1
            if e_3 == i:
                n = 2
            if n != -1:
                jkle[j] = [ikle[i][0], ipoin+s_1, ipoin+i]
                jkle[j+1] = [ipoin+s_1, ikle[i][1], ipoin+i]
                j += 2
                neighbours[i][0] = -1
        # ~~> New boundary segments
        if s_2 < 0:
            jkle[j] = [ikle[i][1], ipoin+ielem+k, ipoin+i]
            jkle[3*ielem+k] = [ipoin+i, ipoin+ielem+k, ikle[i][2]]
            meshj[k] = [ikle[i][1], ikle[i][2]]
            j += 1
            k += 1
        # ~~> New inside segment, swapped
        else:
            # look for the connection with s_1
            e_1, e_2, e_3 = neighbours[s_2]
            n = -1
            if e_1 == i:
                n = 0
            if e_2 == i:
                n = 1
            if e_3 == i:
                n = 2
            if n != -1:
                jkle[j] = [ikle[i][1], ipoin+s_2, ipoin+i]
                jkle[j+1] = [ipoin+s_2, ikle[i][2], ipoin+i]
                j += 2
                neighbours[i][1] = -1
        # ~~> New boundary segments
        if s_3 < 0:
            jkle[j] = [ikle[i][2], ipoin+ielem+k, ipoin+i]
            jkle[3*ielem+k] = [ipoin+i, ipoin+ielem+k, ikle[i][0]]
            meshj[k] = [ikle[i][2], ikle[i][0]]
            j += 1
            k += 1
        # ~~> New inside segment, swapped
        else:
            # look for the connection with s_1
            e_1, e_2, e_3 = neighbours[s_3]
            n = -1
            if e_1 == i:
                n = 0
            if e_2 == i:
                n = 1
            if e_3 == i:
                n = 2
            if n != -1:
                jkle[j] = [ikle[i][2], ipoin+s_3, ipoin+i]
                jkle[j+1] = [ipoin+s_3, ikle[i][0], ipoin+i]
                j += 2
                neighbours[i][2] = -1

    # ~~> Building the new mesh and the new neighbouring
    meshx = np.resize(meshx, ipoin+ielem+isbnd)
    meshy = np.resize(meshy, ipoin+ielem+isbnd)
    meshx[ipoin:ipoin+ielem] = np.sum(meshx[meshk], axis=1)/3.
    meshy[ipoin:ipoin+ielem] = np.sum(meshy[meshk], axis=1)/3.
    meshx[ipoin+ielem:] = np.sum(meshx[meshj], axis=1)/2.
    meshy[ipoin+ielem:] = np.sum(meshy[meshj], axis=1)/2.
    neighbours = Triangulation(meshx, meshy, jkle).get_cpp_triangulation()\
                                                  .get_neighbors()
    # ~~> Reset ipobo while you are at it
    jpobo = np.zeros(ipoin+ielem+isbnd, np.int)

    for n in range(ipoin+ielem+isbnd):
        s_1, s_2, s_3 = neighbours[n]
        e_1, e_2, e_3 = jkle[n]
        if s_1 < 0:
            jpobo[e_1] = e_1+1
            jpobo[e_2] = e_2+1
        if s_2 < 0:
            jpobo[e_2] = e_2+1
            jpobo[e_3] = e_3+1
        if s_3 < 0:
            jpobo[e_3] = e_3+1
            jpobo[e_1] = e_1+1

    return jkle, meshx, meshy, jpobo, meshj, meshk

def trace_ray2xy(ikle, meshx, meshy, neighbours, e_i, xyi, e_n, xyn):
    """
    This assumes that you cannot go back on your ray.
    """
    # ~~> latest addition to the ray
    a_x, b_x, c_x = meshx[ikle[e_n]]
    a_y, b_y, c_y = meshy[ikle[e_n]]
    b_i = get_barycentric_weights(xyi, (a_x, a_y), (b_x, b_y), (c_x, c_y))
    pnt = {'n':1,
           'xy':[xyi],
           'e':[e_n],
           'b':[b_i],
           'd':[np.power(xyi[0]-xyn[0], 2) + np.power(xyi[1]-xyn[1], 2)]}

    # ~~> convergence on distance to target xyn
    accuracy = np.power(10.0,
                        -5+np.floor(np.log10(abs(a_x+b_x+c_x+a_y+b_y+c_y))))
    if pnt['d'][0] < accuracy:
        return True, pnt

    # ~~> get the ray through to the farthest neighbouring edges
    k_s = []
    d_s = []
    for k in [0, 1, 2]:
        xyj = get_segment_intersection(\
                (meshx[ikle[e_n][k]], meshy[ikle[e_n][k]]),
                (meshx[ikle[e_n][(k+1)%3]],
                 meshy[ikle[e_n][(k+1)%3]]),
                xyi, xyn)
        if xyj == []:
            # there are n_o intersection with that edges
            continue
        e_j = neighbours[e_n][k]
        if e_j == e_i:
            # you should not back track on your ray
            continue
        xyj = xyj[0]
        dij = np.power(xyi[0]-xyj[0], 2) + np.power(xyi[1]-xyj[1], 2)
        k_s.append(k)
        d_s.append(dij)
    if d_s != []:
        k = k_s[np.argmax(d_s)]
        e_j = neighbours[e_n][k]
        xyj = get_segment_intersection(\
                (meshx[ikle[e_n][k]], meshy[ikle[e_n][k]]),
                (meshx[ikle[e_n][(k+1)%3]],
                 meshy[ikle[e_n][(k+1)%3]]),
                xyi, xyn)[0]
        djn = np.power(xyn[0]-xyj[0], 2) + np.power(xyn[1]-xyj[1], 2)

        # ~~> Possible recursive call
        if True or djn > accuracy:    # /!\ this may be a problem
            if e_j < 0:
                # you have reach the end of the line
                b_j = get_barycentric_weights(xyj, (a_x, a_y),
                                              (b_x, b_y), (c_x, c_y))
                pnt['n'] += 1
                pnt['xy'].insert(0, xyj)
                pnt['e'].insert(0, e_n)
                pnt['b'].insert(0, b_j)
                pnt['d'].insert(0, djn)
                return djn < accuracy, pnt
            else:
                found, ray = trace_ray2xy(ikle, meshx, meshy, neighbours, e_n,
                                          xyj, e_j, xyn)
                ray['n'] += 1
                ray['xy'].append(xyi)
                ray['e'].append(e_n)
                ray['b'].append(b_i)
                ray['d'].append(dij)
                return found, ray

    # ~~> convergence on having found the appropriate triangle
    b_n = is_inside_triangle(xyn, (a_x, a_y), (b_x, b_y), (c_x, c_y))
    if b_n != []:
        pnt['n'] += 1
        pnt['xy'].insert(0, xyn)
        pnt['e'].insert(0, e_n)
        pnt['b'].insert(0, b_n)
        pnt['d'].insert(0, 0.0)
        return True, pnt

    # ~~> you should not be here !
    return False, pnt

def xys_locate_mesh(xyo, ikle, meshx, meshy, tree=None, neighbours=None):

    # ~~> get to the nearest element
    oet = -1
    obr = [0.0, 0.0, 0.0]
    e_o, b_o, tree = near_locate_mesh(np.array(xyo), ikle, meshx, meshy, tree)
    if b_o == []:
        found, ray, neighbours = xy_trace_mesh(\
                e_o, [np.sum(meshx[ikle[e_o]])/3.0,
                      np.sum(meshy[ikle[e_o]])/3.0],
                xyo, ikle, meshx, meshy, neighbours)
        if found:
            obr = ray['b'][ray['n']]
            oet = ray['e'][ray['n']]
    else:
        obr = b_o
        oet = e_o

    if oet == -1:
        return [-1, -1, -1], obr
    return ikle[oet], obr

def cross_mesh(polyline, ikle, meshx, meshy, tree=None, neighbours=None):
    # ~~ Intersection nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ipt = []
    iet = []
    ibr = []

    # ~~ Locate nodes of the polyline ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xyo = polyline[0]
    e_o, b_o, tree = near_locate_mesh(xyo, ikle, meshx, meshy, tree)

    for i in range(len(polyline)-1):
        xyi = polyline[i+1]
        e_i, b_i, tree = near_locate_mesh(xyi, ikle, meshx, meshy, tree)
        if b_o == [] and b_i == []:
            rank = 0
            found, e_j, xyj, _ = dicho_locate_mesh(rank, e_o, xyo, e_i, xyi,
                                                   ikle, meshx, meshy, tree)
            if not found:
                raise TelemacException(\
                   '... Could not find easily a_n intersection with the mesh')
            found, rayo, neighbours = xy_trace_mesh(e_j, xyj, xyo,
                                                    ikle, meshx, meshy,
                                                    neighbours)
            for j in range(rayo['n'])[:-1]:
                ipt.append(rayo['xy'][j])
                iet.append(rayo['e'][j])
                ibr.append(rayo['b'][j])
            found, rayi, neighbours = xy_trace_mesh(e_j, xyj, xyi,
                                                    ikle, meshx, meshy,
                                                    neighbours)
            for j in range(rayi['n'])[(rayi['n']-1)::-1]:
                ipt.append(rayi['xy'][j])
                iet.append(rayi['e'][j])
                ibr.append(rayi['b'][j])
        elif b_i == [] and b_o != []:
            found, rayi, neighbours = xy_trace_mesh(e_o, xyo, xyi,
                                                    ikle, meshx, meshy,
                                                    neighbours)
            for j in range(rayi['n'])[(rayi['n']-1)::-1]:
                ipt.append(rayi['xy'][j])
                iet.append(rayi['e'][j])
                ibr.append(rayi['b'][j])
        elif b_i != [] and b_o == []:
        # it is necessary to reverse the ray for a case with first end outside
            found, rayo, neighbours = xy_trace_mesh(e_i, xyi, xyo,
                                                    ikle, meshx, meshy,
                                                    neighbours)
            for j in range(rayo['n']): #[(rayo['n']-1)::-1]:
                ipt.append(rayo['xy'][j])
                iet.append(rayo['e'][j])
                ibr.append(rayo['b'][j])
        else:
            found, rayi, neighbours = xy_trace_mesh(e_o, xyo, xyi,
                                                    ikle, meshx, meshy,
                                                    neighbours)
            for j in range(rayi['n'])[(rayi['n']-1)::-1]:
                ipt.append(rayi['xy'][j])
                iet.append(rayi['e'][j])
                ibr.append(rayi['b'][j])

        xyo = xyi
        b_o = b_i
        e_o = e_i

    return (ipt, iet, ibr), tree, neighbours


def slice_mesh(polyline, ikle, meshx, meshy, tree=None):
    """
    A new method to slice through a triangular mesh (replaces cross_mesh)
    """
    xys = []
    douplets = []
    # ~~> Calculate the minimum mesh resolution
    dxy = math.sqrt(min(np.square(np.sum(np.fabs(\
                          meshx[ikle]-meshx[np.roll(ikle, 1)]), axis=1)/3.0) + \
                         np.square(np.sum(np.fabs(\
                          meshy[ikle]-meshy[np.roll(ikle, 1)]), axis=1)/3.0)))
    accuracy = np.power(10.0, -8+np.floor(np.log10(dxy)))

    xyo = np.array(polyline[0])
    for i in range(len(polyline)-1):
        xyi = np.array(polyline[i+1])
        dio = math.sqrt(sum(np.square(xyo-xyi)))

        # ~~> Resample the line to that minimum mesh resolution
        rsmpline = np.dstack((np.linspace(xyo[0], xyi[0], num=int(dio/dxy)),
                              np.linspace(xyo[1], xyi[1], num=int(dio/dxy))))[0]
        nbpoints = len(rsmpline)
        nbneighs = min(8, len(ikle))
        # ~~> Filter closest 8 elements (please create a good mesh) a_s a halo
        # around the polyline
        halo = np.zeros((nbpoints, nbneighs), dtype=np.int)
        for i in range(nbpoints):
            _, elem = tree.query(rsmpline[i], nbneighs)
            halo[i] = elem
        halo = np.unique(halo)

        # ~~> Get the intersecting halo (on a smaller mesh connectivity)
        edges = Triangulation(meshx, meshy, ikle[halo]).get_cpp_triangulation()\
                                                       .get_edges()

        # ~~> Last filter, all nodes that are on the polyline
        olah = []
        nodes = np.unique(edges)
        for node in nodes:  # TODO(jcp): replace b_y numpy calcs
            if get_distance_point_to_line((meshx[node],
                                           meshy[node]), xyo, xyi) < accuracy:
                olah.append(node)
        ijsect = list(zip(olah, olah))
        xysect = [(meshx[i], meshy[i]) for i in olah]
        lmsect = [(1.0, 0.0) for i in range(len(ijsect))]
        mask = np.zeros((len(edges), 2), dtype=bool)
        for i in olah:
            mask = np.logical_or(edges == i, mask)
        edges = np.compress(np.logical_not(np.any(mask, axis=1)), edges, axis=0)

        # ~~> Intersection with remaining edges
        for edge in edges:
            xyj = get_segment_intersection((meshx[edge[0]], meshy[edge[0]]),
                                           (meshx[edge[1]], meshy[edge[1]]),
                                           xyo, xyi)
            if xyj != []:
                ijsect.append(edge)     # nodes from the mesh
                xysect.append(tuple(xyj[0]))   # intersection (xo,yo)
                lmsect.append((xyj[1], 1.0-xyj[1]))   # weight along each each

        # ~~> Final sorting along keys x and y
        xysect = np.array(xysect, dtype=[('x', '<f4'), ('y', '<f4')])
        xysort = np.argsort(xysect, order=('x', 'y'))

        # ~~> Move on to next point
        for i in xysort:
            xys.append(xysect[i])
            douplets.append((ijsect[i], lmsect[i]))
        xyo = xyi

    return xys, douplets

def tessellate_poly(i2s, debug=True):
    """
    A new method to tessellate a series of polygons, the biggest one
        being the outter domain, while the smaller ones are islands.
    The method start b_y building a Delaunay triangulation of all points
        before correcting and swapping those edges that d_o not belong.
    """

    # ~~ Ensuring Clockwisines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if debug:
        print('     +> Ensuring clockwisiness of input polygons')
    if debug:
        print('        - Distinguishing the outter polygon')
    pol_areas = np.zeros(len(i2s.poly), dtype=np.float)
    if debug:
        pbar = ProgressBar(maxval=i2s.npoin).start()
    for ipoly in range(len(i2s.poly)):
        if i2s.type[ipoly] == 1:
            x, y = i2s.poly[ipoly].T
            pol_areas[ipoly] = abs(np.sum(x[:-1]*y[1:] - x[1:]*y[:-1]))/2.0
        if debug:
            pbar.update(ipoly)
    if debug:
        pbar.finish()
    pol_main = np.argmax(pol_areas)
    if debug:
        print('        - Forcing clockwise / anti-clockwise prerties')
    if debug:
        pbar = ProgressBar(maxval=i2s.npoin).start()
    if i2s.type[pol_main] == 1:
        if is_clockwise(i2s.poly[pol_main]):
            i2s.poly[pol_main] = np.flipud(i2s.poly[pol_main])
        for ipoly in range(len(i2s.poly)):
            if i2s.type[ipoly] == 1 and ipoly != pol_main:
                if not is_clockwise(i2s.poly[ipoly]):
                    i2s.poly[ipoly] = np.flipud(i2s.poly[ipoly])
            if debug:
                pbar.update(ipoly)
    if debug:
        pbar.finish()

    # ~~ Sorting out polygons' edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if debug:
        print('     +> Sorting out egdes of the input polygons')
    meshx = np.zeros(i2s.npoin)
    meshy = np.zeros(i2s.npoin)
    pol_edges = []
    pol_nbres = np.zeros(i2s.npoin, dtype=np.int)
    cpt = 0
    if debug:
        pbar = ProgressBar(maxval=i2s.npoin).start()
    values = np.arange(len(i2s.poly), dtype=np.int)[np.argsort(pol_areas)][::-1]
    for ipoly in values:
        edges = []
        # number the edges
        lpt = len(i2s.poly[ipoly])
        meshx[cpt:cpt+lpt] = i2s.poly[ipoly].T[0][:lpt]
        meshy[cpt:cpt+lpt] = i2s.poly[ipoly].T[1][:lpt]
        pol_nbres[cpt:cpt+lpt] = ipoly
        if i2s.type[ipoly] == 0:
            edges = np.arange(len(i2s.poly[ipoly]), dtype=np.int) + cpt
        elif i2s.type[ipoly] == 1:
            edges = np.arange(len(i2s.poly[ipoly])+1, dtype=np.int) + cpt
            edges[-1] = edges[0]
        cpt += len(i2s.poly[ipoly])
        if debug:
            pbar.update(cpt)
        pol_edges.extend(zip(edges[:-1], edges[1:]))
        pol_areas[ipoly] = 0.0
        if i2s.type[ipoly] == 1:
            pol_areas[ipoly] = \
            abs(np.sum(meshx[edges[:-1]]*meshy[edges[1:]] - \
                          meshx[edges[1:]]*meshy[edges[:-1]]))/2.0
    if debug:
        pbar.finish()
    pol_main = np.argmax(pol_areas)
    pol_edges = np.array(zip(*np.sort(np.asarray(pol_edges)).T),
                         dtype=[('h', int), ('t', int)])
    pol_edges = pol_edges[np.argsort(pol_edges, order=('h', 't'))]

    # ~~ Delaunay triangulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if debug:
        print('     +> Delaunay triangulation of all points '\
              '(outside and inside)')
    ikle2 = Delaunay(np.dstack((meshx, meshy))[0]).vertices

    # ~~ Constraining the triangulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if debug:
        print('     +> Constraining the triangulation to the input polygon')
    # Note:
    #   Because the nodes in the polygon are ordered, it should be easy to
    #   identify when a_n ( edges[0] != edges[0]-1 ), and check whether it is
    #   a poly-closing edge or a_n edge to be swapped.
    #   /!\ that get_edges() is ordered, but not necessarily oriented

    if debug:
        print('        - Extracting edges from the triangulation')
    tri = Triangulation(meshx, meshy, ikle2).get_cpp_triangulation()
    tri_neigh = tri.get_neighbors()
    tri_edges = np.asarray(zip(*np.sort(tri.get_edges()).T),
                           dtype=[('h', int), ('t', int)])
    tri_edges = tri_edges[np.argsort(tri_edges, order=('h', 't'))]

    if debug:
        print('        - Searching for missing edges in the triangulation')
    ipe = 0
    ite = 0
    swp_edges = []
    if debug:
        pbar = ProgressBar(maxval=len(pol_edges)).start()
    while ipe < len(pol_edges):
        edge = pol_edges[ipe]
        k = np.searchsorted(tri_edges['h'][ite:], edge[0])
        ite += k
        if debug:
            pbar.update(ipe)
        if ite == len(tri_edges):     # remaining edges not in the triangulation
            swp_edges.extend(pol_edges[ipe:])
            ipe = len(pol_edges)
        elif edge == tri_edges[ite]:  # found it first trial
            ite += 1
            ipe += 1
        elif edge[0] < tri_edges[ite][0]:
            swp_edges.append(edge)
            ipe += 1
            ite = np.searchsorted(tri_edges['h'], edge[0])
        else:
            ite += 1
    if debug:
        pbar.finish()

    if debug:
        print('        - Swapping the edges to follow input polygons')
    # Note:
    #   A ray is traced along the missing edges s_o to break through and rebuild
    #   the triangulation incorporating the missing edge. If there is only one
    #   edge crossing the missing edge, then a simple swap is done. Otherwise,
    #   it re-create two separate triangulations on either side of the edge
    #   n_1 and n_2 are the two nodes of the input polygon defining the missing
    #   edge in the triangulation
    #   e_i and e_n are the elements including n_1 and n_2 respectively, the
    #   opposit segement of which intersect with the missing edge
    cpt = 0
    if debug:
        pbar = ProgressBar(maxval=len(swp_edges)).start()
    for n_1, n_2 in swp_edges:

        eis = np.extract(np.any(np.equal(n_1, ikle2), axis=1),
                         np.arange(len(ikle2)))
        for e_i in eis:
            ki1, ki2 = np.setxor1d([n_1], ikle2[e_i])
            seg1 = mplPath.Path([[meshx[n_1], meshy[n_1]],
                                 [meshx[n_2], meshy[n_2]]])
            seg2 = mplPath.Path([[meshx[ki1], meshy[ki1]],
                                 [meshx[ki2], meshy[ki2]]])
            if seg1.intersects_path(seg2) == 1:
                break
        ens = np.extract(np.any(np.equal(n_2, ikle2), axis=1),
                         np.arange(len(ikle2)))
        for e_n in ens:
            kn1, kn2 = np.setxor1d([n_2], ikle2[e_n])
            seg1 = mplPath.Path([[meshx[n_1], meshy[n_1]],
                                 [meshx[n_2], meshy[n_2]]])
            seg2 = mplPath.Path([[meshx[kn1], meshy[kn1]],
                                 [meshx[kn2], meshy[kn2]]])
            if seg1.intersects_path(seg2) == 1:
                break
        if debug:
            pbar.update(cpt)
        cpt += 1
        # ~~> if n_o segment crossing, continue (the problem may have corrected
        # itself)
        if len(np.intersect1d(eis, ens) == 2):
            continue
        # ~~> if only one segment crossing, swap it
        elif e_i in tri_neigh[e_n] and e_n in tri_neigh[e_i]:
            if is_ccw((meshx[n_1], meshy[n_1]),
                      (meshx[n_2], meshy[n_2]),
                      (meshx[ki1], meshy[ki1])):
                ikle2[e_i] = [n_1, n_2, ki1]
            else:
                ikle2[e_i] = [ki1, n_2, n_1]
            if is_ccw((meshx[n_1], meshy[n_1]),
                      (meshx[n_2], meshy[n_2]),
                      (meshx[ki2], meshy[ki2])):
                ikle2[e_n] = [n_1, n_2, ki2]
            else:
                ikle2[e_n] = [ki2, n_2, n_1]
            ejs = [e_i, e_n]
        # ~~> if more than one segment crossing, rebuild two local
        # triangulations on either side of the missing edge
        else:
            ejs = [e_i]     # elements sliced b_y the missing edge
            ki1s = [n_1, ki1]   # nodes on one side
            ki2s = [n_1, ki2]   # nodes on the other side
            e_j = e_i
            while e_j != e_n:
                # the missing edge goes between ki1 and ki2 shared between e_i
                # and e_j
                for e_j in tri_neigh[e_i]:
                    if e_j < 0:
                        continue
                    if ki1 in ikle2[e_j] and ki2 in ikle2[e_j]:
                        break
                ejs.append(e_j)
                # there shall be only one
                kj0 = np.setxor1d([ki1, ki2], ikle2[e_j])[0]
                seg1 = mplPath.Path([[meshx[n_1], meshy[n_1]],
                                     [meshx[n_2], meshy[n_2]]])
                seg2 = mplPath.Path([[meshx[ki1], meshy[ki1]],
                                     [meshx[kj0], meshy[kj0]]])
                if seg1.intersects_path(seg2) == 1:
                    if kj0 in ki2s:
                        raise TelemacException(\
                                'mesh.p_y : tessellate_poly function')
                    ki2s.append(kj0)
                    ki2 = kj0
                else:
                    if kj0 in ki1s:
                        raise TelemacException(\
                                'mesh.p_y : tessellate_poly function')
                    ki1s.append(kj0)
                    ki1 = kj0
                e_i = e_j
            if n_2 not in ki1s:
                ki1s.append(n_2)
            ki1s = np.array(ki1s)
            if n_2 not in ki2s:
                ki2s.append(n_2)
            ki2s = np.array(ki2s)

            # ~~> new triangulations on either side of the missing edge
            # side one: assumed convex
            sub1_ikle2 = Delaunay(np.dstack((meshx[ki1s],
                                             meshy[ki1s]))[0]).vertices
            # if side one not convex
            isoxy = np.column_stack((np.sum(meshx[ki1s[sub1_ikle2]],
                                            axis=1)/3.0,
                                     np.sum(meshy[ki1s[sub1_ikle2]],
                                            axis=1)/3.0))
            poly_path = mplPath.Path(np.dstack((meshx[ki1s], meshy[ki1s]))[0])
            sub1_ikle2 = sub1_ikle2[poly_path.contains_points(isoxy)]
            # side two: assumed convex
            sub2_ikle2 = Delaunay(np.dstack((meshx[ki2s],
                                             meshy[ki2s]))[0]).vertices
            # if side two not convex
            isoxy = np.column_stack((np.sum(meshx[ki2s[sub2_ikle2]],
                                            axis=1)/3.0,
                                     np.sum(meshy[ki2s[sub2_ikle2]],
                                            axis=1)/3.0))
            poly_path = mplPath.Path(np.dstack((meshx[ki2s], meshy[ki2s]))[0])
            sub2_ikle2 = sub2_ikle2[poly_path.contains_points(isoxy)]

            # ~~> like for like replacement
            ikle2[ejs] = np.concatenate((ki1s[sub1_ikle2], ki2s[sub2_ikle2]))

        # ~~> finding out the new neighbouring state
        halo = np.unique(np.concatenate((ejs, tri_neigh[ejs].ravel())))
        insiders = {}
        bounders = {}
        for elem, i in zip(ikle2[halo], halo):
            n_k = bounders.keys()
            for k in [0, 1, 2]:
                if (elem[k], elem[(k+1)%3]) not in n_k:
                    bounders.update({(elem[(k+1)%3], elem[k]):i})
                else:
                    j = bounders[(elem[k], elem[(k+1)%3])]
                    insiders.update({(elem[k], elem[(k+1)%3]):[i, j]})
                    del bounders[(elem[k], elem[(k+1)%3])]
        for elem, i in zip(ikle2[halo], halo):
            for k in [0, 1, 2]:
                if (elem[k], elem[(k+1)%3]) in insiders:
                    e_a, e_b = insiders[(elem[k], elem[(k+1)%3])]
                    if e_a == i:
                        tri_neigh[i][k] = e_b
                    if e_b == i:
                        tri_neigh[i][k] = e_a
                if (elem[(k+1)%3], elem[k]) in insiders:
                    e_a, e_b = insiders[(elem[(k+1)%3], elem[k])]
                    if e_a == i:
                        tri_neigh[i][k] = e_b
                    if e_b == i:
                        tri_neigh[i][k] = e_a

    if debug:
        pbar.finish()

    # ~~ Removing outside elements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if debug:
        print('     +> Removing the outside elements')

    if debug:
        print('        - Computing iso-barycentres')
    isoxy = np.column_stack((np.sum(meshx[ikle2], axis=1)/3.0,
                             np.sum(meshy[ikle2], axis=1)/3.0))

    if debug:
        print('        - Identifying and removing who is out')
    # pol_main is the domain boundary. The others are islands
    tri_inside = np.ones(len(ikle2), dtype=bool)
    tri_mask = np.arange(len(ikle2), dtype=np.int)\
              [np.equal(*pol_nbres[np.sort(ikle2).T[0::2]])]
    pkle2 = pol_nbres[ikle2.T[0]]
    cpt = 0
    if debug:
        pbar = ProgressBar(maxval=len(i2s.poly)).start()
    values = np.arange(len(i2s.poly), dtype=np.int)[np.argsort(pol_areas)][::-1]
    for ipoly in values:
        poly_path = mplPath.Path(i2s.poly[ipoly])
        pol_mask = np.equal(ipoly, pkle2[tri_mask])
        tri = tri_mask[pol_mask]
        if ipoly == pol_main:
            tri_inside[tri] = poly_path.contains_points(isoxy[tri])
        else:
            tri_inside[tri] = np.logical_not(\
                    poly_path.contains_points(isoxy[tri]))
        tri_mask = tri_mask[np.logical_not(pol_mask)]
        if debug:
            pbar.update(cpt)
        cpt += 1
    if debug:
        pbar.finish()
    ikle2 = ikle2[tri_inside]

    if debug:
        print('     +> set ipobo')
    ipob2 = np.arange(i2s.npoin, dtype=np.int) + 1

    return ikle2, ipob2, meshx, meshy

def filter_mesh_resolution(meshx, meshy, ikle2, resolut, factor, debug=True):

    tree = cKDTree(np.column_stack(meshx, meshy))

    if debug:
        pbar = ProgressBar(maxval=len(resolut)).start()
    while True:
        found = 0
        for i_o in range(len(resolut)):
            xyo = [meshx[i_o], meshy[i_o]]
            for elem, i in zip(*tree.query(xyo, 50)):
                if resolut[i] > resolut[i_o] + elem * factor:
                    found += 1
                    resolut[i] = resolut[i_o] + elem * factor
            if debug:
                pbar.update(max(len(resolut)-found, 0))
        if found == 0:
            break
    if debug:
        pbar.finish()

def filter_poly_resolution(i2s, factor, resolut, debug=True):

    _, _, meshx, meshy = tessellate_poly(i2s, debug)
    tree = cKDTree(np.column_stack(meshx, meshy))

    if debug:
        pbar = ProgressBar(maxval=len(i2s.npoin)).start()
    while True:
        found = 0
        for i_o in range(len(i2s.vals)):
            xyo = [meshx[i_o], meshy[i_o]]
            for elem, i in zip(*tree.query(xyo, 50)):
                if resolut[i] > resolut[i_o] + elem * factor:
                    found += 1
                    resolut[i] = resolut[i_o] + elem * factor
            if debug:
                pbar.update(max(len(resolut)-found, 0))
        if found == 0:
            break
    if debug:
        pbar.finish()


def cross_check_boundaries(meshx, meshy, ikle2, ipob0, debug=True):

    if debug:
        print('     +> removing cross boundaries')
    # ~~> establish neighborhood
    _ = Triangulation(meshx, meshy, ikle2).get_cpp_triangulation()\
                                                  .get_neighbors()
    crosspts = np.where(ipob0 == 4)[0]
    sout = - np.ones(len(crosspts), dtype=np.int)
    eout = np.extract(np.any(ipob0[ikle2] == 4, axis=1), np.arange(len(ikle2)))
    for i in range(len(crosspts)):
        if debug:
            print('        - node: ', crosspts[i], ' ( ',
                  meshx[crosspts[i]], ' , ', meshy[crosspts[i]], ' )')
        jout = ikle2[eout[np.where(ikle2[eout] == crosspts[i])[0]]]
        if len(jout) == 2:
            kout = 0
        else: # how about if you have more than 3 nodes connecting ?
            kout = 0
            if len(np.intersect1d(jout[0], jout[1])) == 2:
                if len(np.intersect1d(jout[1], jout[2])) == 1:
                    kout = 2
            elif len(np.intersect1d(jout[1], jout[2])) == 1:
                kout = 1
        sout[i] = eout[np.where(ikle2[eout] == crosspts[i])[0][kout]]
        # removing only one of the elements
    mask = np.ones(len(ikle2), dtype=np.bool)
    mask[sout] = False

    return ikle2[mask]

def get_ipobo(meshx, meshy, ikle2, debug=True):

    if debug:
        print('     +> sorting boundary nodes')
    # ~~> establish neighborhood
    tri_neigh = Triangulation(meshx, meshy, ikle2).get_cpp_triangulation()\
                                                  .get_neighbors()
    ebounds = []
    if debug:
        print('        - identify')
    pbar = ProgressBar(maxval=len(ikle2)).start()
    for i in range(len(ikle2)):
        if tri_neigh[i, 0] < 0:
            ebounds.append([ikle2[i][0], ikle2[i][1]])
        if tri_neigh[i, 1] < 0:
            ebounds.append([ikle2[i][1], ikle2[i][2]])
        if tri_neigh[i, 2] < 0:
            ebounds.append([ikle2[i][2], ikle2[i][0]])
        pbar.update(i)
    pbar.finish()
    # ~~> assemble the enssemble of boundary segments
    if debug:
        print('        - assemble')
    pbounds = join_segments(ebounds)
    # ~~> define ipobo from a_n arbitrary start point
    if debug:
        print('        - set')
    ipob2 = np.zeros(len(meshx), dtype=np.int)
    iptfr = 0
    for point in pbounds:
        for n in point[1:]:
            iptfr += 1
            ipob2[n] = iptfr

    return ipob2, pbounds

def show_boundary_nodes(meshx, meshy, ikle2, debug=True):

    if debug:
        print('     +> highlight boundary nodes')
    # ~~> establish neighborhood
    tri_neigh = Triangulation(meshx, meshy, ikle2).get_cpp_triangulation()\
                                                  .get_neighbors()
    # ~~> check boundary nodes
    ipob0 = np.zeros(len(meshx), dtype=np.int)
    pbar = ProgressBar(maxval=len(ikle2)).start()
    for i in range(len(ikle2)):
        if tri_neigh[i, 0] < 0:
            ipob0[ikle2[i][0]] += 1
            ipob0[ikle2[i][1]] += 1
        if tri_neigh[i, 1] < 0:
            ipob0[ikle2[i][1]] += 1
            ipob0[ikle2[i][2]] += 1
        if tri_neigh[i, 2] < 0:
            ipob0[ikle2[i][2]] += 1
            ipob0[ikle2[i][0]] += 1
        pbar.update(i)
    pbar.finish()

    return ipob0

def show_node_connections(meshx, meshy, ikle2, debug=True):

    if debug:
        print('     +> count node connections')
    # ~~> establish neighborhood
    tri_edges = Triangulation(meshx, meshy, ikle2).get_cpp_triangulation()\
                                                  .get_edges()
    # ~~> count node connections
    ipob1 = np.zeros(len(meshx), dtype=np.int)
    for edge in tri_edges:
        ipob1[edge] += 1

    return ipob1

def merge_min_4_nodes(meshx, meshy, ikle2, where, debug=True):

    if debug:
        print('     +> removing nodes connected to 4 or less others')
    # ~~> preparing masks
    mask_ikle2 = np.ones(len(ikle2), dtype=bool)
    mask_npoin = np.ones(len(meshx), dtype=bool)
    for ipoin in where:
        # ~~> elements surounding the problematic node
        w = np.where(ikle2 == ipoin)[0]
        # ~~> nodes surounding the problematic node
        n_1 = np.setdiff1d(np.unique(ikle2[w]), [ipoin])[0]
        # ~~> two elements are to be masked
        for elem in w:
            if n_1 in ikle2[elem]:
                mask_ikle2[elem] = False
            else:
                for i in range(3):
                    if ikle2[elem][i] == ipoin:
                        ikle2[elem][i] = n_1
        mask_npoin[ipoin] = False # you need renumbering
        # TODO: search for the opposing edges before joining -- or swap if
        # needed

    return mask_npoin, ikle2[mask_ikle2]


def cleave_max_7_nodes(meshx, meshy, ikle2, where, debug=True):
    # TODO: finich implementation ?

    raise TelemacException('Implementation not finished')
    if debug:
        print('     +> cleaving nodes connected of ')
    # ~~> preparing masks
    add_ikle2 = []
    add_meshx = []
    add_meshy = []
    for ipoin in where:
        # ~~> elements surounding the problematic node
        w = np.where(ikle2 == ipoin)[0]
        edges = Triangulation(meshx, meshy, ikle2[w]).get_cpp_triangulation()\
                                                     .get_edges()
        edge = edges[0]
        # ~~> average distance between all nodes
        dist = 0.
        d_x = 0.
        d_y = 0.
        for edge in edges:
            dxmin = max(meshx[edge]) - min(meshx[edge])
            dymin = max(meshy[edge]) - min(meshy[edge])
            dist += math.sqrt(dxmin*dxmin + dymin*dymin)/len(edges)
            d_x += dxmin/len(edges)
            d_y += dymin/len(edges)
        print(repr((dist, d_x, d_y)))
        # ~~> nodes surounding the problematic node
        n = np.setdiff1d(np.unique(ikle2[w]), [ipoin])
        x_c = np.sum(meshx[n])/len(n)
        y_c = np.sum(meshy[n])/len(n)
        print(repr((x_c, y_c, x_c+dist*d_x/2./math.sqrt(d_x*d_x+d_y*d_y),
                    y_c+dist*d_y/2./math.sqrt(d_x*d_x+d_y*d_y))))


        dxmin = min(meshx[n][1:] - meshx[n][:-1])
        dymin = min(meshy[n][1:] - meshy[n][:-1])
        dist = math.sqrt(dxmin*dxmin + dymin*dymin)
        dxmax = abs(max(meshx[n][1:] - meshx[n][:-1]))
        dymax = abs(max(meshy[n][1:] - meshy[n][:-1]))
        print(repr((ipoin, x_c, y_c, dist, dxmax, dymax)))
        # ~~> two elements are to be added
        #print( ipoin,w,n dd)
    return None, None

def remove_duplicate_nodes(meshx, meshy, ikle2, alpha, debug=True):
    # /!\ just checking ...
    # TODO: implement the removal / merge of the duplicated nodes

    if debug:
        print('     +> checking for duplicated nodes')
    tree = cKDTree(np.column_stack((meshx, meshy)))

    if debug:
        pbar = ProgressBar(maxval=len(meshx)).start()
    for j in range(len(meshx)): # TODO: parallel implementation

        # ~~> subset of points (one closest nodes should be OK)
        i = np.setdiff1d(tree.query([meshx[j], meshy[j]], 2)[1], [j])[0]
        dist = math.sqrt((meshx[j]-meshx[i])*(meshx[j]-meshx[i]) + \
                              (meshy[j]-meshy[i])*(meshy[j]-meshy[i]))
        if dist < alpha:
            print('        - distance between nodes:' + str(j) + \
                    ' and ' + str(i) + ' is ' + repr(dist), j)

        if debug:
            pbar.update(j)
    if debug:
        pbar.finish()

    return meshx, meshy, ikle2


def remove_extra_nodes(meshx, meshy, ikle2, debug=True):

    if debug:
        print('     +> removing extra nodes')
    # ~~> knolg(npoin3) gives the global node number such that
    #   for i = 1,npoin3: Fwrite(i) = Fread(knolg(i)) and is ordered
    knolg = np.sort(np.unique(np.ravel(ikle2), return_index=True)[0])
    knogl = dict(zip(knolg, range(len(knolg))))
    likle = - np.ones_like(ikle2, dtype=np.int)
    if debug:
        print('        - renumbering connectivity')
        pbar = ProgressBar(maxval=len(ikle2)).start()
    for k in range(len(ikle2)):
        likle[k] = [knogl[ikle2[k][0]],
                    knogl[ikle2[k][1]],
                    knogl[ikle2[k][2]]]
        if debug:
            pbar.update(k)
    if debug:
        pbar.finish()
        print('        - removing extra nodes from mesh')
        pbar = ProgressBar(maxval=len(ikle2)).start()
    for k in knolg:
        meshx[knogl[k]] = meshx[k]
        meshy[knogl[k]] = meshy[k]
        if debug:
            pbar.update(k)
    if debug:
        pbar.finish()

    return meshx[:len(knolg)], meshy[:len(knolg)], likle

def map_thin_plate_spline(meshx, meshy, bathx, bathy, bathz, npoin,
                          alpha, debug=True):

    if debug:
        print('     +> creating a proximity tree')
    tree = cKDTree(np.column_stack((bathx, bathy)))
    # ~~> mask for the proximity points
    mask = np.zeros(npoin, dtype=np.int)

    varsor = np.zeros(len(meshx), dtype=np.float)

    if debug:
        print('     +> thin plate spline mapping')
        pbar = ProgressBar(maxval=len(meshx)).start()
    mat_a = np.zeros((npoin+3, npoin+3), dtype=np.float)
    row = np.zeros(npoin+3, dtype=np.float)
    for j in range(len(meshx)): # TODO: parallel implementation
        # ~~> subset of points
        mask = tree.query([meshx[j], meshy[j]], npoin)[1]   # indices
        # ~~> weight
        dt = 0.
        for i in range(npoin):
            dt += np.sum(\
               2*np.sqrt(np.power(bathx[mask][i]-bathx[mask][i+1:], 2)+\
                         np.power(bathy[mask][i]-bathy[mask][i+1:], 2)))
        dt = dt / (npoin*npoin)

        # ~~> main matrix
        mat_a = mat_a * 0.              # resetting mat_a
        # ~~> upper triangular elements
        for i in range(npoin):
            row = row * 0.       # reseting row
            x_i = bathx[mask][i] * np.ones(npoin-(i+1))
            y_i = bathy[mask][i] * np.ones(npoin-(i+1))
            r_a = np.sqrt(np.square(bathx[mask][(i+1):]-x_i)+\
                         np.square(bathy[mask][(i+1):]-y_i))
            r_a[np.where(r_a == 0)[0]] = 0.001
            rowslice = np.square(r_a)*np.log(r_a) # / 2.0 ?
            row[i] = alpha * dt * dt
            row[(i+1):(i+1+len(rowslice))] = rowslice
            row[(npoin):] = np.array([1, bathx[mask][i], bathy[mask][i]])
            mat_a[i] = row

        # ~~> build vec_a symetric matrix of mat_a
        mat = np.mat(mat_a)
        mat = mat + mat.T

        # ~~> second member
        vec_b = np.append(bathz[mask], [0, 0, 0])

        # ~~> linear algebra
        vec_a = np.linalg.solve(mat, vec_b)

        # ~~> resulting interpolation
        r_a = np.sqrt(np.square(bathx[mask]-meshx[j]*np.ones(npoin))+\
                         np.square(bathy[mask]-meshy[j]*np.ones(npoin)))
        # /!\ must check for NaN from calculating log(0) ...
        # Thank you David Roscoe !
        varsor[j] = vec_a[-3]+vec_a[-2]*meshx[j]+vec_a[-1]*meshy[j]+\
                        np.nansum(np.square(r_a)*vec_a[:npoin]*np.log(r_a))

        if debug:
            pbar.update(j)
    if debug:
        pbar.finish()

    return varsor


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban"
__date__ = "$12-Dec-2012 08:51:29$"

def main():

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Command line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+'~'*72+'\n')
    parser = ArgumentParser(\
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
Tools for editing, sampling and interpolating through triangular meshes.
        '''),
        usage=' (--help for help)\n---------\n      =>  '\
                '%(prog)s [options] in-file.slf out-file.slf\n---------')
    parser.add_argument(\
        "-x", "--crosses", action="store_true", dest="xpts", default=False,
        help="check and remove cross boundary points")
    parser.add_argument(\
        "-4", "--minnode", action="store_true", dest="mnod", default=False,
        help="check and merge interior nodes with less than 5 neighbours")
    parser.add_argument(\
        "-8", "--maxnode", action="store_true", dest="xnod", default=False,
        help="check and cleave interior nodes with more than 7 neighbours")
    parser.add_argument(\
        "--dupnode", action="store_true", dest="dupnode", default=False,
        help="check and remove duplicated nodes")
    parser.add_argument(\
        "--map", dest="mapping", default=None,
        help="name of the points for mapping (thin plate spline)")
    parser.add_argument(\
        "-n", "--cloudsize", dest="npoin", default=None,
        help="size of the proximity cloud (thin plate spline)")
    parser.add_argument(\
        "-r", "--smoothing", dest="alpha", default=None,
        help="smoothing factor (thin plate spline)")
    parser.add_argument(\
        "args", metavar='2x selafin files', default='', nargs=2,
        help="first file: your input mesh; second file: the processed mesh")
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads code name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    slf_file = options.args[0]
    if not path.exists(slf_file):
        raise TelemacException(\
            '... the provided (first) file does not seem to exist: {} \n\n'\
            .format(slf_file))

    out_file = options.args[1]
    slf = Selafin(slf_file)

    # ~~> Monitoring changes
    updt_ipobo = False
    updt_npoin = False

    # ~~> check and remove boundary nodes where boundary lines intersect
    if options.xpts:
        # ipob0 is 0 for inside node, or the noumber of boundary segment
        # attached
        # np.sum((ipob0[slf.ikle2]>0),axis=1) > 1 highlights elements with at
        # least two boundary nodes
        ipob0 = show_boundary_nodes(slf.meshx, slf.meshy, slf.ikle2, debug=True)
        if max(ipob0) > 2:
            slf.ikle2 = cross_check_boundaries(slf.meshx, slf.meshy, slf.ikle2,
                                               ipob0, debug=True)
            updt_ipobo = True
            slf.ikle3 = slf.ikle2
            slf.nelem2 = len(slf.ikle2)
            slf.nelem3 = slf.nelem2

    # ~~> remove inside nodes that have 4 or less connecting nodes
    if options.mnod:
        # ipob0 is 0 for inside node, or the noumber of boundary segment
        # attached
        ipob0 = show_boundary_nodes(slf.meshx, slf.meshy, slf.ikle2, debug=True)
        # ipob1 counts the number of connections around a node
        ipob1 = show_node_connections(slf.meshx, slf.meshy,
                                      slf.ikle2, debug=True)
        if np.any(np.logical_and(ipob1 < 5, ipob0 == 0)):
            _, slf.ikle2 = merge_min_4_nodes(\
                      slf.meshx, slf.meshy, slf.ikle2,
                      np.where(np.logical_and(ipob1 < 5, ipob0 == 0))[0],
                      debug=True)
            updt_ipobo = True
            updt_npoin = True
            slf.ikle3 = slf.ikle2
            slf.nelem2 = len(slf.ikle2)
            slf.nelem3 = slf.nelem2

    # ~~> cleave inside nodes that have 8 or more connecting nodes
    if options.xnod:
        # ipob0 is 0 for inside node, or the noumber of boundary segment
        # attached
        ipob0 = show_boundary_nodes(slf.meshx, slf.meshy, slf.ikle2, debug=True)
        # ipob1 counts the number of connections around a node
        ipob1 = show_node_connections(slf.meshx, slf.meshy,
                                      slf.ikle2, debug=True)
        if np.any(np.logical_and(ipob1 > 7, ipob0 == 0)):
            _, slf.ikle2 = cleave_max_7_nodes(\
                      slf.meshx, slf.meshy, slf.ikle2,
                      np.where(np.logical_and(ipob1 > 7, ipob0 == 0))[0],
                      debug=True)
            updt_ipobo = True
            updt_npoin = True
            slf.ikle3 = slf.ikle2
            slf.nelem2 = len(slf.ikle2)
            slf.nelem3 = slf.nelem2

    # ~~> check and remove duplicated nodes
    if options.dupnode:
        alpha = 1.0 # minimum resolution allowed
        if options.alpha:
            alpha = float(options.alpha)
        slf.meshx, slf.meshy, slf.ikle2 = remove_duplicate_nodes(\
                  slf.meshx, slf.meshy, slf.ikle2, alpha, debug=True)
        updt_ipobo = True
        updt_npoin = True
        slf.npoin2 = len(slf.meshx)
        slf.npoin3 = slf.npoin2
        slf.ikle3 = slf.ikle2

    # ~~> remove singled out nodes, and renumber ikle accordingly
    if updt_npoin:
        slf.meshx, slf.meshy, slf.ikle2 = remove_extra_nodes(\
                  slf.meshx, slf.meshy, slf.ikle2, debug=True)
        slf.npoin2 = len(slf.meshx)
        slf.npoin3 = slf.npoin2
        slf.ikle3 = slf.ikle2

    if updt_ipobo:
        # ~~> Effective way to limit the search and sort for boundary nodes
        # ipob0 is 0 form inside node, or the noumber of boundary segment
        # attached
        # np.sum((ipob0[slf.ikle2]>0),axis=1) > 1 highlights elements with at
        # least two boundary nodes
        ipob0 = show_boundary_nodes(slf.meshx, slf.meshy, slf.ikle2, debug=True)
        mask = np.sum((ipob0[slf.ikle2] > 0), axis=1) > 1
        slf.ipob2, _ = get_ipobo(slf.meshx, slf.meshy, slf.ikle2[mask],
                                 debug=True)
        slf.ipob3 = slf.ipob2

    varsor = np.zeros(len(slf.meshx), dtype=np.float)
    if options.mappping:
        mapp = Selafin(options.mapping)
        npoin = min(20, mapp.npoin2)
        if options.npoin:
            npoin = min(int(options.npoin), mapp.npoin2)
        alpha = 0.2
        if options.alpha:
            alpha = float(options.alpha)
        varsor = map_thin_plate_spline(\
                  slf.meshx, slf.meshy, mapp.meshx, mapp.meshy,
                  mapp.get_variables_at(0, [0])[0], npoin, alpha, debug=True)

    # ~~> new Selafin writer
    slf.fole = {}
    slf.fole.update({'hook':open(out_file, 'wb')})
    slf.fole.update({'name':out_file})
    slf.fole.update({'endian':">"})     # big endian
    slf.fole.update({'float':('f', 4)})  # single precision

    print('     +> Write Selafin header')
    slf.append_header_slf()

    print('     +> Write Selafin core')
    slf.append_core_time_slf(0.0)
    slf.append_core_vars_slf([varsor])
    slf.fole['hook'].close()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)

if __name__ == "__main__":
    main()
