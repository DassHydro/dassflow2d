#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   @brief Function for plot.py
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
import argparse
import matplotlib.pyplot as plt
# Unsued but must be kept
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
# ~~> dependencies towards standard python
from postel.plot_vnv import vnv_plot2d, vnv_plot1d_polylines, \
        vnv_plot1d_history
from postel.plot2d import plot2d_spectrum
from postel.plot1d import plot1d
from utils.exceptions import TelemacException

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#


def add_options_fig(parser):
    """
    Add options for a figure (x_label, y_label, title, fig_name...)

    @param parser (ArgumentParser) The parser

    @returns the updated parser
    """
    # Options to save the file instead of displaying it
    parser.add_argument(
        "-f", "--figure-name",
        dest="fig_name", default="",
        help="If given the figure will be saved in fig_name instead "
             "of beeing displayed")

    return parser


def add_options_var(subparser):
    """
    Defines options for var action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'var', help='Plot a scalar map for a given variable and time/record')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="Name of the variable to display")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=0,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")
    parser.add_argument(
        "--mesh", action="store_true",
        dest="mesh", default=False,
        help="Adds the mesh to the display")
    parser = add_options_fig(parser)

    return subparser


def arg_points(string):
    """
    Definition of a point for argparse
    """
    n_coords = string.count(",") + 1

    if n_coords == 2:
        try:
            x, y = map(float, string.split(','))
            return x, y
        except Exception:
            raise argparse.ArgumentTypeError("Points must be x,y")
    elif n_coords == 3:
        try:
            x, y, z = map(float, string.split(','))
            return x, y, z
        except Exception:
            raise argparse.ArgumentTypeError("Points must be x,y,z")
    else:
        raise argparse.ArgumentTypeError("Points must be either x,y or x,y,z")


def add_options_timeseries_on_polyline(subparser):
    """
    Defines options for time_poly action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('time_poly',
                                  help='Plot a timeseries on a polyline')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="Name of the variable to display")
    parser.add_argument(
        "--poly",
        dest="poly", type=arg_points, nargs='+',
        help="List of points for the polyline x,y space separated")

    parser = add_options_fig(parser)

    return subparser


def add_options_history(subparser):
    """
    Defines options for history action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('history',
                                  help='Plot history for a list of points')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="Name of the variable to display")
    parser.add_argument(
        "--points",
        dest="points", type=arg_points, nargs='+',
        help="List of points for which to plot history x,y space separated")

    parser = add_options_fig(parser)

    return subparser


def add_options_mesh2d(subparser):
    """
    Defines options for mesh2 action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'mesh2d',
        help='Plot a 2d mesh can add boundary or liquid boundary info')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    # the boundary file option
    parser.add_argument(
        "-b", "--boundary-file",
        dest="bnd_file", default=None,
        help="Name of the boundary file")

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--bnd", action="store_true",
        dest="bnd", default=False,
        help="Adding type of boundary for each boundary node")
    group.add_argument(
        "--liq-bnd", action="store_true",
        dest="liq_bnd", default=False,
        help="Adding number of liquid boundary for each boundary node")

    parser = add_options_fig(parser)

    return subparser


def add_options_3d_scalar_map(subparser):
    """
    Defines options for 3d_scalar_map action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        '3d_scalar_map',
        help='Plot a 3d representation of a 2D variable '
        'using variable values as Z coordinates')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")

    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="Name of the variable to display")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=-1,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")
    parser.add_argument(
        "--mesh", action="store_true",
        dest="mesh", default=False,
        help="Adds the mesh to the display")
    parser = add_options_fig(parser)

    return subparser


def list_of_points(arg):
    """
    Change the string argument in a 2 dimension array.

    @param arg (string)

    @returns a list of list
    """
    res = [[float(v) for v in r.lstrip('(').rstrip(')').split(',')]
           for r in arg.replace(' ', '').split(';')]
    return res


def add_options_vertical_slice(subparser):
    """
    Defines options for v_slice action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'v_slice',
        help='Plot a vertical slice of a 3d mesh along a polyline')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")

    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="Name of the variable to display")

    parser.add_argument("--poly",
                        dest='poly',
                        type=arg_points,
                        nargs='+',
                        help="Choose the  points (xi,yi) where to extract \
                              use --poly='x1,y1 x2,y2...'")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=-1,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")
    parser.add_argument(
        "--mesh", action="store_true",
        dest="mesh", default=False,
        help="Adds the mesh to the display")
    parser = add_options_fig(parser)
    return subparser


def add_options_horizontal_slice(subparser):
    """
    Defines options for h_slice action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'h_slice',
        help='Plot a horizontal slice of a 3d mesh along a plane number')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")

    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="Name of the variable to display")

    parser.add_argument("--plane",
                        dest='plane',
                        type=int,
                        help="Plane number to slice")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=-1,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")
    parser.add_argument(
        "--mesh", action="store_true",
        dest="mesh", default=False,
        help="Adds the mesh to the display")
    parser = add_options_fig(parser)
    return subparser


def add_options_spe(subparser):
    """
    Defines options for spec action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'spec',
        help='Plot the spectrum of a given point over quandrangle mesh')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    parser.add_argument(
        "-p", "--point",
        dest="point", type=int,
        help="Number of the point to display")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=-1,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")

    parser = add_options_fig(parser)

    return subparser


def add_options_spe_freq(subparser):
    """
    Defines options for spec-freq action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'spec-freq',
        help='Plot the frequency of the spectrum of given points')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    parser.add_argument(
        "-p", "--points",
        dest="points", default=None, type=int, nargs='*',
        help="Number of the points to display if none given display them all")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=-1,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")

    parser = add_options_fig(parser)

    return subparser


def add_options_spe_ang(subparser):
    """
    Defines options for spec-ang action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser(
        'spec-ang',
        help='Plot the angular dispersion of the spectrum of given points')

    parser.add_argument(
        "input_file", default=None,
        help="Name of the input file extension also defines the input format")
    parser.add_argument(
        "-p", "--points",
        dest="points", default=None, type=int, nargs='*',
        help="Number of the points to display if none given display them all")
    parser.add_argument(
        "-r", "--record",
        dest="record", type=int, default=-1,
        help="Record to display (If -1 is given will return the last one)")
    parser.add_argument(
        "-t", "--time",
        dest="time", type=float, default=None,
        help="Time to display (will take the closest record)")

    parser = add_options_fig(parser)

    return subparser


def add_options_poly(subparser):
    """
    Defines options for poly action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('poly',
                                  help='Plot history for a list of points')

    parser.add_argument(
        "shp_file",
        help="Name of the shape file")

    return subparser



def plot_var(res, var, record=-1, time=None, add_mesh=False,
             fig_name=''):
    """
    Plot a scalar map for the given variable and time record

    @param res (TelemacFile) Structure to file from which data will be read
    @param var (str) Name of the variable to plot
    @param record (str) Record to plot
    @param time (str) If >= 0.0 will get nearest record to that time (This
    overwrites record)
    @param add_mesh (boolean) If True overlay the mesh on top of the scalar map
    @param fig_name (str) If not empty save the plot in that file instead of
    showing it
    """
    # If time is positive searched for record
    if time is not None:
        record = res.get_closest_record(time)
        time = res.times[record]
    else:
        time = res.times[record]

    if var not in res.varnames:
        raise TelemacException("{} is not in :\n{}".format(var, res.varnames))

    vnv_plot2d(var, res, plot_mesh=add_mesh, record=record,
               filled_contours=True,
               aspect_ratio="equal",
               fig_name=fig_name)


def plot_mesh2d(res, display_bnd=False,
                display_liq_bnd=False, fig_name=''):
    """
    Plot a 2d triangular mesh with either boundary conditions or liquid
    boundary number

    @param res (TelemacFile) File from wich to read the mesh
    @param display_bnd (boolean) If True display boundary type for each
        boundary node
    @param display_liq_bnd (boolean) If True display liquidi boundary number
        for each boundary node
    @param fig_name (str) If not empty save the plot in that file instead of
    showing it
    """
    if (display_bnd or display_liq_bnd) and res.boundary_file == '':
        raise TelemacException(
            "bnd_file is mandatory if using --bnd or --liq-bnd")

    try:
        varname = res.varnames[0]
    except IndexError:
        varname = ''
    vnv_plot2d(varname, res,
               plot_mesh=True,
               annotate_bnd=display_bnd,
               annotate_liq_bnd=display_liq_bnd,
               aspect_ratio="equal",
               fig_name=fig_name)


def plot_timeseries_on_polyline(res, var_name, poly, records=None,
                                fig_name=''):
    """
    Plot a value over a polyline for a range of time

    @param res (TelemacFile) Struct of the file from which to extract the data
    @param var_name (str) Name of the variable fro which to extract the data
    @param poly (list) List of polyline points
    @param records (list) List of record for which to extrac the data
    @param fig_name (str) If not empty saving in that file
    """

    if records is None:
        records = range(res.ntimestep)

    vnv_plot1d_polylines(var_name, res, poly, record=records,
                         fig_name=fig_name)


def plot_history(res, var_name, points, fig_name=''):
    """
    Plot values of points over a range of records

    @param res (TelemacFile) Struct of the file from which to extract the data
    @param var_name (str) Name of the variable fro which to extract the data
    @param points (list) List of points (x,y) in 2d (x,y,z) in 3d
    @param fig_name (str) If not empty saving in that file
    """
    vnv_plot1d_history(var_name, res, points=points, fig_name=fig_name)


def plot_vertical_slice(res, varname, poly,
                        record=-1, time=None, add_mesh=False,
                        fig_name=''):
    """
    Plot a vertical slice of a 3d mesh

    @param res (TelemacFile) Struct to file from which data will be read
    @param varname (str) Name of the variable to plot
    @param poly (list) List of polyline points
    @param record (str) Record to plot
    @param time (str) If >= 0.0 will get nearest record to that time (This
    overwrites record)
    @param add_mesh (boolean) If True overlay the mesh on top of the scalar map
    @param fig_name (str) If not empty save the plot in that file instead of
    showing it
    """
    # If time is positive searched for record
    if time is not None:
        rrecord = res.get_closest_record(time)
    else:
        rrecord = record

    vnv_plot2d(
        varname,
        res,
        poly=poly,
        plot_mesh=add_mesh,
        record=rrecord,
        filled_contours=True,
        aspect_ratio="equal",
        fig_name=fig_name)


def plot_horizontal_slice(res, varname, plane, record=-1, time=None,
                          add_mesh=False, fig_name=''):
    """
    Plat an horizontal slice of a 3d mesh for given plane number

    @param res (TelemacFile) Struct to file from which data will be read
    @param varname (str) Name of the variable to plot
    @param plane (int) Number of the plane from which to extrac
    @param record (str) Record to plot
    @param time (str) If >= 0.0 will get nearest record to that time (This
    overwrites record)
    @param add_mesh (boolean) If True overlay the mesh on top of the scalar map
    @param fig_name (str) If not empty save the plot in that file instead of
    showing it
    """
    # If time is positive searched for record
    if time is not None:
        rrecord = res.get_closest_record(time)
    else:
        rrecord = record

    if plane < 0:
        iplane = res.nplan + plane
    else:
        iplane = plane

    vnv_plot2d(varname, res, plane=iplane, plot_mesh=add_mesh,
               filled_contours=True,
               aspect_ratio="equal",
               record=rrecord, fig_name=fig_name)


def plot_spe(res, point, record=-1, time=None, fig_name=''):
    """
    Plotting a specter for a given point
    """
    # If time is given searched for record
    if time is not None:
        record = res.get_closest_record(time)
        time = res.times[record]
    else:
        time = res.times[record]

    # Getting name of the variable containing the spectrum for point 2
    var_name = res.get_spectrum_varname(point)

    # Getting value of the spectrum
    data = res.get_data_value(var_name, record)

    fig, axe = plt.subplots(figsize=(12, 10))

    # Ploting it
    plot2d_spectrum(fig, axe, res.meshx, res.meshy, res.ikle2, data)

    axe.set_title("At time {}".format(time))

    if fig_name != '':
        print(" "*8+"~> Plotting {}".format(fig_name))
        fig.savefig(fig_name)
    else:
        plt.show()

    plt.close('all')


def plot_spe_freq(res, points=None, record=-1, time=None, fig_name=''):
    """
    Plotting the frequency for the spectrum on points
    """
    # If time is given searched for record
    if time is not None:
        record = res.get_closest_record(time)
        time = res.times[record]
    else:
        time = res.times[record]

    if points is None:
        points = res.get_list_spectrum_points()

    fig, axe = plt.subplots()

    for point in points:
        # Getting list of frequencies and spectrum value
        freq, spectrum = res.get_spectrum(point, record)
        # Plotting it
        plot1d(axe, freq, spectrum, plot_label='point {:06d}'.format(point),
               x_label='Frequencies', y_label='Spectrum')

    axe.legend()

    axe.set_title("At time {}".format(time))

    if fig_name != '':
        print(" "*8+"~> Plotting {}".format(fig_name))
        fig.savefig(fig_name)
    else:
        plt.show()

    plt.close('all')


def plot_spe_ang(res, points=None, record=-1, time=None, fig_name=''):
    """
    Plotting the angular dispersion for the spectrum on points
    """
    # If time is given searched for record
    if time is not None:
        record = res.get_closest_record(time)
        time = res.times[record]
    else:
        time = res.times[record]

    if points is None:
        points = res.get_list_spectrum_points()

    fig, axe = plt.subplots(subplot_kw={'projection': 'polar'})

    for point in points:
        # Getting list of frequencies and spectrum value
        theta, disp = res.get_angular_dispersion(point, record, radian=True)
        # Plotting it
        plot1d(axe, theta, disp, plot_label='point {:06d}'.format(point),
               x_label='Angles', y_label='Angular dispersion')
        axe.set_theta_zero_location("N")
        axe.set_theta_direction(-1)

    axe.legend()

    axe.set_title("At time {}".format(time))

    if fig_name != '':
        print(" "*8+"~> Plotting {}".format(fig_name))
        fig.savefig(fig_name)
    else:
        plt.show()

    plt.close('all')

def plot_poly(polys):
    """
    Plotting polygons
    """
    # my_poly is a list of list containing the (x,y) coordinate
    if np.array(polys).ndim == 2:
        ppolys = [polys]
    else:
        ppolys = polys
    print('Number of shapes: %d' % len(polys))

    # First polyline
    _, ax = plt.subplots()
    for poly in ppolys:
        ppoly = np.array(poly)
        ax.plot(ppoly[:, 0], ppoly[:, 1], 'o-')

    ax.grid()

    plt.show()
