#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   @brief check_mesh tools
"""
from argparse import ArgumentParser
from math import degrees, acos
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from os import path, environ

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.tri import TriAnalyzer

from data_manip.extraction.telemac_file import TelemacFile


def check_non_ortho_for_dual_mesh(mesh, plot_hist=False, plot_mesh=False):
    """

    Check non-orthogality of edge on the dual mesh

    @param mesh  (TelemacFile) mesh structure
    @param plot_hist  (Bool) plot histograms?
    @param plot_mesh  (Bool) 2D plot results?

    @returns (double)
    """
    # TODO Don't Work
    print("Mesh quality: \n \
           Check non-orthogality of edge on the dual mesh")

    ikle2 = mesh.ikle2
    npoin2 = mesh.npoin2
    nelem2 = mesh.nelem2
    # fill edges and neighbours
    mesh.set_mpl_tri()

    nedges = np.shape(mesh.edges)[0]
    ifacel = np.zeros((nedges, 2), dtype=np.int64)
    cfve = np.zeros(2, dtype=np.float64)
    cocg = np.zeros((nelem2, 2), dtype=np.float64)
    c_node_1 = np.zeros(2, dtype=np.float64)
    c_node_2 = np.zeros(2, dtype=np.float64)
    # f_ortho = np.zeros(nedges, dtype=np.float64)
    elt_ortho = np.zeros(nelem2, dtype=np.float64) + 2
    node_ortho = np.zeros(npoin2, dtype=np.float64) + 2

    for i_elt in range(nelem2):
        node_1 = ikle2[i_elt, 0]
        node_2 = ikle2[i_elt, 1]
        node_3 = ikle2[i_elt, 2]
        cocg[i_elt, 0] = (mesh.meshx[node_1]
                          + mesh.meshx[node_2]
                          + mesh.meshx[node_3]) / 3.
        cocg[i_elt, 1] = (mesh.meshy[node_1]
                          + mesh.meshy[node_2]
                          + mesh.meshy[node_3]) / 3.
    ifac = 0
    for i_elt in range(nelem2):
        neig_1 = mesh.neighbours[i_elt, 0]
        neig_2 = mesh.neighbours[i_elt, 1]
        neig_3 = mesh.neighbours[i_elt, 2]
        if i_elt < neig_1 or neig_1 == -1:
            ifacel[ifac, 0] = i_elt
            ifacel[ifac, 1] = neig_1
            ifac += 1

        if i_elt < neig_2 or neig_2 == -1:
            ifacel[ifac, 0] = i_elt
            ifacel[ifac, 1] = neig_2
            ifac += 1
        if i_elt < neig_3 or neig_3 == -1:
            ifacel[ifac, 0] = i_elt
            ifacel[ifac, 1] = neig_3
            ifac += 1

    for i_edge in range(nedges):
        iel_1 = ifacel[i_edge, 0]
        iel_2 = ifacel[i_edge, 1]
        node_1 = mesh.edges[i_edge, 0]
        node_2 = mesh.edges[i_edge, 1]
        c_node_1[0] = mesh.meshx[node_1]
        c_node_1[1] = mesh.meshy[node_1]
        c_node_2[0] = mesh.meshx[node_2]
        c_node_2[1] = mesh.meshy[node_2]

        if iel_1 >= 0 and iel_2 >= 0:
            cfcg = 0.5 * (c_node_1 + c_node_2)

            cfve[0] = -(c_node_2[1] - c_node_1[1])
            cfve[1] = c_node_2[0] - c_node_1[0]

            cfcgve1 = cfcg - cocg[iel_1, :]
            cfcgve2 = cfcg - cocg[iel_2, :]

            ccgve = cocg[iel_2, :] - cocg[iel_1, :]

            p_norm_orthof1 = np.linalg.norm(cfve) * np.linalg.norm(cfcgve1)
            term1_orthof1 = abs(np.vdot(cfve, cfcgve1) / p_norm_orthof1)

            p_norm_orthof2 = np.linalg.norm(cfve) * np.linalg.norm(cfcgve2)
            term1_orthof2 = abs(np.vdot(cfve, cfcgve2) / p_norm_orthof2)

            p_norm_orthof = np.linalg.norm(cfve) * np.linalg.norm(ccgve)
            term2_orthof = abs(np.vdot(cfve, ccgve) / p_norm_orthof)

            # print(i_edge, term1_orthof1, term2_orthof, cfve / np.linalg.norm(cfve))
            f_ortho1 = min(term1_orthof1, term2_orthof)
            f_ortho2 = min(term1_orthof2, term2_orthof)

            # if iel_2==8918 or iel_1 ==8918:
            #     print(i_edge, term1_orthof1, term2_orthof, cfve / np.linalg.norm(cfve))
            elt_ortho[iel_1] = min(elt_ortho[iel_1], f_ortho1)
            elt_ortho[iel_2] = min(elt_ortho[iel_2], f_ortho2)

            node_ortho[node_1] = min(node_ortho[node_1], f_ortho1)
            node_ortho[node_2] = min(node_ortho[node_2], f_ortho2)

    # print(elt_ortho[8918])

    # for i, val in enumerate(elt_ortho):
    #     if val >0.8:
    #         print(i, val)

    s_node_ortho = pd.Series(node_ortho)
    print(s_node_ortho.describe())

    s_elt_ortho = pd.Series(elt_ortho)
    print(s_elt_ortho.describe())

    if plot_hist:
        fig, ax0 = plt.subplots(figsize=(10, 5))
        ax0.set_xlabel("Orthogonality (-)")
        ax0.set_ylabel("Occurences")
        s_node_ortho.plot.hist(bins=5, alpha=1.0)
        plt.show()

    if plot_mesh:
        fig, ax0 = plt.subplots()
        ax0.set_aspect('equal')
        elt_ortho_plot = ax0.tripcolor(mesh.meshx,
                                       mesh.meshy,
                                       mesh.ikle2,
                                       facecolors=elt_ortho,
                                       edgecolors='k')
        fig.colorbar(elt_ortho_plot)
        ax0.set_title('Orthogonality mesh quality criterion')
        ax0.set_xlabel('X (m)')
        ax0.set_ylabel('Y (m)')
        plt.show()

    return node_ortho


def check_circle_ratios_num_of_neig(mesh, plot_hist=False, plot_mesh=False):
    """
    Check the ratios of circle lengh around an elements in a 2D mesh
    Check the number of elements around one node in a 2D mesh

    @param mesh  (TelemacFile) mesh structure
    @param plot_hist  (Bool) plot histograms?
    @param plot_mesh  (Bool) 2D plot results?

    @returns (double) ratios, number of elements around one node
    """

    print("Mesh quality: \n \
           Check the ratios of circle lengh around an element \n \
           Check the number of elements around one node in a 2D mesh")

    ikle2 = mesh.ikle2
    npoin2 = mesh.npoin2
    nelem2 = mesh.nelem2
    tri = mesh.tri
    mesh_circle_ratio = np.zeros(npoin2)
    mesh_circle_w = np.zeros(npoin2)
    c_ratio = TriAnalyzer(tri).circle_ratios()
    s_c_ratio = pd.Series(c_ratio)
    print(s_c_ratio.describe())

    if plot_hist:
        fig, ax0 = plt.subplots(figsize=(10, 5))
        ax0.set_xlabel("Circle ratio (0.5 best) (-)")
        ax0.set_ylabel("Occurences")
        s_c_ratio.plot.hist(bins=5, alpha=1.0)
        plt.draw()

    for i_elt in range(nelem2):
        node_1 = ikle2[i_elt, 0]
        mesh_circle_w[node_1] += 1
        mesh_circle_ratio[node_1] = mesh_circle_ratio[node_1] + c_ratio[i_elt]
        node_2 = ikle2[i_elt, 1]
        mesh_circle_w[node_2] += 1
        mesh_circle_ratio[node_2] = mesh_circle_ratio[node_2] + c_ratio[i_elt]
        node_3 = ikle2[i_elt, 2]
        mesh_circle_w[node_3] += 1
        mesh_circle_ratio[node_3] = mesh_circle_ratio[node_3] + c_ratio[i_elt]

    mesh_circle_ratio /= mesh_circle_w
    s_circle_w = pd.Series(mesh_circle_w)
    print(s_circle_w.describe())

    if plot_hist:
        fig, ax0 = plt.subplots(figsize=(10, 5))
        ax0.set_xlabel("Number of triangles around a node (-)")
        ax0.set_ylabel("Occurences")
        s_circle_w.plot.hist(bins=5, alpha=1.0)
        plt.show()

    if plot_mesh:
        fig, ax0 = plt.subplots()
        ax0.set_aspect('equal')
        c_ratio_plot = ax0.tripcolor(mesh.meshx,
                                     mesh.meshy,
                                     mesh.ikle2,
                                     facecolors=c_ratio,
                                     edgecolors='k')
        fig.colorbar(c_ratio_plot)
        ax0.set_title('Circle ratio mesh quality criterion')
        ax0.set_xlabel('X (m)')
        ax0.set_ylabel('Y (m)')
        plt.show()

    return mesh_circle_ratio, mesh_circle_w


def check_angle_tri(mesh, min_angle=5):
    """
    Check the minimum angle in mesh

    @param mesh  (TelemacFile) mesh structure
    @param min_angle (float) minumum angle

    @returns
    """

    ikle2 = mesh.ikle2
    npoin2 = mesh.npoin2
    nelem2 = mesh.nelem2
    mesh_a = np.zeros(npoin2)
    nb_min_angl_points = 0

    for i_elt in range(nelem2):
        angl_n1 = 0
        angl_n2 = 0
        angl_n3 = 0
        node_1 = ikle2[i_elt, 0]
        node_2 = ikle2[i_elt, 1]
        node_3 = ikle2[i_elt, 2]
        a_len = ((mesh.meshx[node_2] - mesh.meshx[node_3]) ** 2 +
                 (mesh.meshy[node_2] - mesh.meshy[node_3]) ** 2) ** 0.5

        b_len = ((mesh.meshx[node_1] - mesh.meshx[node_3]) ** 2 +
                 (mesh.meshy[node_1] - mesh.meshy[node_3]) ** 2) ** 0.5

        c_len = ((mesh.meshx[node_1] - mesh.meshx[node_2]) ** 2 +
                 (mesh.meshy[node_1] - mesh.meshy[node_2]) ** 2) ** 0.5
        if (a_len + b_len >= c_len) and (b_len + c_len >= a_len) and (c_len + a_len >= b_len):
            angl_n1 = degrees(acos((b_len ** 2 +
                                    c_len ** 2 -
                                    a_len ** 2) / (2 * b_len * c_len)))
            angl_n2 = degrees(acos((c_len ** 2 +
                                    a_len ** 2 -
                                    b_len ** 2) / (2 * c_len * a_len)))
            angl_n3 = degrees(acos((a_len ** 2 +
                                    b_len ** 2 -
                                    c_len ** 2) / (2 * a_len * b_len)))

        if angl_n1 <= min_angle:
            mesh_a[node_1] += 1
            nb_min_angl_points += 1
        if angl_n2 <= min_angle:
            mesh_a[node_2] += 1
            nb_min_angl_points += 1
        if angl_n3 <= min_angle:
            mesh_a[node_3] += 1
            nb_min_angl_points += 1
    if nb_min_angl_points > 0:
        list_node = np.where(mesh_a > 0)[0]
        print("Global number of nodes with angle < {}:".format(min_angle))
        for node in list_node:
            print(node)
    return mesh_a


def check_overconstraint_nodes(mesh):
    """
    Check the number of overcobtraint node in a 2D mesh

    @param mesh  (TelemacFile) mesh structure

    @returns (int) Number of overcontraint nodes
    """

    print("Mesh quality: \n \
           Check number of overconstraint nodes")

    ikle2 = mesh.ikle2
    npoin2 = mesh.npoin2
    nelem2 = mesh.nelem2
    mesh_t = np.zeros(npoin2)
    nb_isolated_points = 0
    nb_overconst_points = 0
    isolated_points = []
    overconst_points = []

    for i_elt in range(nelem2):
        node_1 = ikle2[i_elt, 0]
        mesh_t[node_1] = mesh_t[node_1] + 2.
        node_2 = ikle2[i_elt, 1]
        mesh_t[node_2] = mesh_t[node_2] + 2.
        node_3 = ikle2[i_elt, 2]
        mesh_t[node_3] = mesh_t[node_3] + 2.

    for i_node in range(npoin2):
        if mesh_t[i_node] < 2:
            mesh_t[i_node] = -1
            nb_isolated_points += 1
            isolated_points.append(i_node)
        elif mesh_t[i_node] < 3:
            mesh_t[i_node] = 1
            nb_overconst_points += 1
            overconst_points.append(i_node)
        else:
            mesh_t[i_node] = 0

    if nb_overconst_points > 0:
        print("WARNING: number of overconstraint nodes: {}" \
              .format(nb_overconst_points))
        print("Global number of overcontraint nodes:")
        for i_node in overconst_points:
            print(i_node + 1)
    else:
        print("There is no overcontraint nodes")

    if nb_isolated_points > 0:
        print("ERROR: number of isolated nodes: {}" \
              .format(nb_isolated_points))
        print("Global number of isolated nodes:")
        for i_node in isolated_points:
            print(i_node + 1)


    return nb_overconst_points, overconst_points, mesh_t


def output_mesh_quality_informations(input_res, file_output_name,
                                     data_names, data_values):
    """
    Create file with check mesh informations

    @param input_res  (TelemacFile) mesh structure
    @param file_output_name  (Char) File name for postprocessing data
    @param data_names  (list) Name of check mesh criteria
    @param data_values  (list) Values of check mesh critera (on nodes)

    """
    res_out = TelemacFile(file_output_name, access='w')
    res_out.read_mesh(input_res)
    res_out._ntimestep = 1
    res_out._times = np.array([0.0])

    res_out.nvar = len(data_names)
    res_out.varnames = data_names
    res_out.varunits = [''] * res_out.nvar
    if res_out._ndim == 3:
        res_out._npoin3 = res_out.npoin2
        res_out._nelem3 = res_out.nelem2
        res_out._ndp3 = res_out.ndp2
        res_out._ipob3 = res_out.ipob2
        res_out._ikle3 = res_out.ikle2
        res_out._ndim = 2
        res_out._nplan = 1
        res_out._meshx = res_out.meshx[0:res_out.npoin2]
        res_out._meshy = res_out.meshy[0:res_out.npoin2]
        res_out._meshz = None
        res_out.typ_elem = 10

    res_out._values = np.ones((1, res_out._nvar, res_out.npoin2),
                              dtype=np.float64)
    for i_var in range(res_out._nvar):
        res_out._values[0, i_var, :] = \
            res_out._values[0, i_var, :] * data_values[i_var].flatten()

    res_out.write()


def main():
    """
    Main program
    @return:
    """
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(
        help='Check command to do', dest='command')

    all_parser = subparsers.add_parser("all")
    circl_ratio_parser = subparsers.add_parser("circl_ratio")
    overconstraint_parser = subparsers.add_parser("overconstraint")
    # orthogonal_parser = subparsers.add_parser("orthogonality")
    check_angle_parser = subparsers.add_parser("check_angle")
    list_parsers = [all_parser,
                    circl_ratio_parser,
                    overconstraint_parser,
                    check_angle_parser]
                    # orthogonal_parser,

    for pars in list_parsers:
        pars.add_argument(
            'file_to_check',
            help="File to check",
            type=str)

        pars.add_argument(
            '--outfile',
            dest="outfile",
            help="outfile name containing the fields of quality mesh",
            type=str,
            default='geo_quality_mesh.slf')

    circl_ratio_parser.add_argument(
        '--plot_hist',
        default=False,
        action="store_true",
        dest='plot_hist',
        help='plot histogramm for circle ratio value')

    circl_ratio_parser.add_argument(
        '--plot_mesh',
        default=False,
        action="store_true",
        dest='plot_mesh',
        help='plot mesh elements for circle ratio field')

    check_angle_parser.add_argument(
        '--angle',
        default=5,
        type=float,
        dest='min_angle',
        help='Minimum angle(deg) to check (5 by default)')

    all_parser.add_argument(
        '--plot_hist_circl_ratio',
        default=False,
        dest='plot_hist_circl_ratio',
        action="store_true",
        help='plot histogramm for circle ratio value')

    all_parser.add_argument(
        '--plot_mesh_circl_ratio',
        default=False,
        action="store_true",
        dest='plot_mesh_circl_ratio',
        help='plot mesh elements for circle ratio field')

    all_parser.add_argument(
        '--angle',
        default=5,
        type=float,
        dest='min_angle',
        help='Minimum angle(deg) to check (5 by default)')

    # orthogonal_parser.add_argument( \
    #         '--plot_hist',
    #         default=False,
    #         dest='plot_hist',
    #         action="store_true",
    #         help='plot histogramm for orthogonality quality element')
    # orthogonal_parser.add_argument( \
    #         '--plot_mesh',
    #         default=False,
    #         action="store_true",
    #         dest='plot_mesh',
    #         help='plot mesh elements for orthogonality quality field')
    # all_parser.add_argument( \
    #         '--plot_hist_ortho',
    #         default=False,
    #         dest='plot_hist',
    #         action="store_true",
    #         help='plot histogramm for orthogonality quality element')
    # all_parser.add_argument( \
    #         '--plot_mesh_ortho',
    #         default=False,
    #         action="store_true",
    #         dest='plot_mesh',
    #         help='plot mesh elements for orthogonality quality field')
    #

    args = parser.parse_args()

    check_dico = {"circl_ratio": False,
                  "check_angle": False,
                  "overconstraint": False,
                  "orthogonal": False}
    if args.command == 'all':
        for key in check_dico:
            check_dico[key] = True
        # no work
        check_dico["orthogonal"] = False
    else:
        check_dico[args.command] = True

    input_file = TelemacFile(args.file_to_check)
    _data_names = []
    _data_values = []

    if check_dico["circl_ratio"]:
        cond_hist = False
        cond_mesh = False
        if args.command == 'all':
            if args.plot_hist_circl_ratio:
                cond_hist = True
            if args.plot_mesh_circl_ratio:
                cond_mesh = True
        else:
            if args.plot_hist:
                cond_hist = True
            if args.plot_mesh:
                cond_mesh = True
        circle_ratio_field, elts_arnd_nodes_field = \
            check_circle_ratios_num_of_neig(input_file,
                                            plot_hist=cond_hist,
                                            plot_mesh=cond_mesh)
        _data_values.append(circle_ratio_field)
        _data_names.append("CIRCLE RATIO")
        _data_values.append(elts_arnd_nodes_field)
        _data_names.append("ELTS ARND NODES")

    if check_dico["check_angle"]:
        mesh_min_angle_field = check_angle_tri(input_file, args.min_angle)
        _data_values.append(mesh_min_angle_field)
        _data_names.append("MIN ANGLE NODES")
    if check_dico["overconstraint"]:
        nb_overconst_nodes, overconst_nodes, over_data_field = \
            check_overconstraint_nodes(input_file)
        _data_values.append(over_data_field)
        _data_names.append("OVER CONST NODES")
    if check_dico["orthogonal"]:
        cond_hist = False
        cond_mesh = False
        if args.command == 'all':
            if args.plot_hist_ortho:
                cond_hist = True
            if args.plot_mesh_ortho:
                cond_mesh = True
        else:
            if args.plot_hist:
                cond_hist = True
            if args.plot_mesh:
                cond_mesh = True
        node_ortho_field = check_non_ortho_for_dual_mesh(input_file,
                                                         plot_hist=cond_hist,
                                                         plot_mesh=cond_mesh)
        _data_values.append(node_ortho_field)
        _data_names.append("NODE ORTHO")
    if _data_names:
        path_dir = path.realpath(args.file_to_check)
        path_outfile = path.join(path.dirname(path_dir), args.outfile)
        output_mesh_quality_informations(input_file,
                                         path_outfile,
                                         _data_names,
                                         _data_values)
    del input_file


def test():
    """
    Test programm
    @return:
    """
    _file_name = path.join(environ['HOMETEL'],
                           'examples',
                           'telemac2d',
                           'gouttedo',
                           'geo_gouttedo.slf')
    # 'telemac2d',
    # 'malpasset',
    # 'geo_malpasset-large.slf')

    _file_output_name = path.join(environ['HOMETEL'],
                                  'examples',
                                  'telemac2d',
                                  'gouttedo',
                                  'geo_gouttedo_out.slf')

    min_angle = 5  # degree
    # Opening the file with TelemacFile
    res = TelemacFile(_file_name)
    node_ortho_field = check_non_ortho_for_dual_mesh(res, False, False)

    circle_ratio_field, elts_arnd_nodes_field = \
        check_circle_ratios_num_of_neig(res)

    mesh_min_angle_field = check_angle_tri(res, min_angle)

    nb_overconst_nodes, overconst_nodes, over_data_field = \
        check_overconstraint_nodes(res)
    _data_names = ["OVER CONST NODES", "CIRCLE RATIO",
                   "ELTS ARND NODES", "NODE ORTHO",
                   "MIN ANGLE NODES"]
    _data_values = [over_data_field, circle_ratio_field,
                    elts_arnd_nodes_field, node_ortho_field,
                    mesh_min_angle_field]
    output_mesh_quality_informations(res,
                                     _file_output_name,
                                     _data_names,
                                     _data_values)

    res.close()


if __name__ == '__main__':
    main()
    # test()
