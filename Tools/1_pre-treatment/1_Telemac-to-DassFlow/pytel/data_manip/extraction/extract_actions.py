"""
Function for extract.py
"""
import argparse
import numpy as np
from utils.exceptions import TelemacException
from data_manip.extraction.telemac_file import TelemacFile


def add_options_csv(parser):
    """
    Add options for a figure (x_label, y_label, title, fig_name...)

    @param parser (ArgumentParser) The parser

    @returns the updated parser
    """
    # Options to save the file instead of displaying it
    parser.add_argument(
        dest="csv_name",
        help="Name of the output csv file")

    parser.add_argument(
        "--delimiter",
        dest="delimiter", default=";",
        help="Delimiter in the csv file")

    return parser


def arg_points(string):
    """
    Definition of a point for argparse
    """
    n_coords = string.count(",") + 1

    if n_coords == 2:
        try:
            x, y = map(float, string.split(','))
            return x, y
        except Exception as e:
            raise argparse.ArgumentTypeError("Points must be x,y")
    elif n_coords == 3:
        try:
            x, y, z = map(float, string.split(','))
            return x, y, z
        except Exception as e:
            raise argparse.ArgumentTypeError("Points must be x,y,z")
    else:
        raise argparse.ArgumentTypeError("Points must be either x,y or x,y,z")


def add_options_timeseries(subparser):
    """
    Defines options for timeseries action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('timeseries',
                                  help='Extract timeseries over\
                                        points or nodes')

    parser.add_argument(dest="file_name",
                        help="Telemac file to extract from")

    parser.add_argument(
        "-v", "--var",
        dest="var", default="",
        help="List of names of the variable to display ',' separated. "
             "If more than one is given csv_name will be appended the "
             "name of the variable")

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--points", default=None,
        dest="points", type=arg_points, nargs='+',
        help="List of points for extraction x,y or x,y,z space separated")

    group.add_argument(
        "--nodes", default=None,
        dest="nodes", type=int, nargs='+',
        help="List of nodes for extraction nodes space separated")

    add_options_csv(parser)

    return subparser


def extract_timeseries(file_name, var_name, points=None, nodes=None):
    """
    Extract timeseries informations on a list of nodes or points

    @param file_name (str) Name of the file from which to extract
    @param var_name (str) Name of the variable to extract
    @param points (List) List of points on which to extract
    @param nodes (List) List of nodes on which to extract

    @returns (List, numpy.array) List of strings (name for each column), the
    extracted data
    """

    res = TelemacFile(file_name)
    header = ['time (s)']

    if points is not None:
        tmp_data = res.get_timeseries_on_points(var_name, points)
        for point in points:
            header.append(str(point))
    elif nodes is not None:
        tmp_data = res.get_timeseries_on_nodes(var_name, nodes)
        for node in nodes:
            header.append(str(node))
    else:
        raise TelemacException("Give at least points or nodes")

    data = np.vstack((res.times, tmp_data))

    res.close()

    return header, data.T


def add_options_mesh2d(subparser):
    """
    Defines options for mesh action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('mesh2d',
                                  help='Extract x y from file')

    parser.add_argument(dest="file_name",
                        help="Telemac file to extract from")

    add_options_csv(parser)

    return subparser


def extract_mesh2d(file_name):
    """
    Extract mesh coordinates

    @param file_name (str) Name of the file from which to extract

    @returns (List, numpy.array) List of strings (name for each column), the
    extracted data
    """

    res = TelemacFile(file_name)
    header = ['X', 'Y']

    data = np.column_stack((res.meshx, res.meshy))

    return header, data


def add_options_spectrum(subparser):
    """
    Defines options for spectrum action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('spectrum',
                                  help='Extract spectrum of a\
                                        given node at a given record/time')

    parser.add_argument(dest="file_name",
                        help="Telemac file to extract from")

    parser.add_argument(
        "-p", "--point",
        dest="point", type=int,
        help="Point of the spectrum to extract")

    parser.add_argument(
        "--radian",
        action='store_true',
        dest="radian", default=False,
        help="If given theta is in radian instead of degree")

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "-r", "--record", default=0,
        dest="record", type=int,
        help="Record of the spectrum to extract")

    group.add_argument(
        "-t", "--time", default=None,
        dest="time", type=float,
        help="Time of the spectrum to extract")

    add_options_csv(parser)

    return subparser


def extract_spectrum(file_name, point, radian=False, record=0, time=None):
    """
    Extract timeseries informations on a list of nodes or points

    @param file_name (str) Name of the file from which to extract
    @param point (int) Point number of the spectrum to extract
    @param radian (bool) If true theta will be given in radian instead
        of degree
    @param record (int) Record to extract
    @param time (float) Time to extract

    @returns (List, numpy.array) List of strings (name for each column), the
    extracted data
    """

    res = TelemacFile(file_name)

    # Getting list of frequencies
    freqs, _ = res.get_spectrum_freq()

    nfreq = len(freqs)
    ntheta = res.npoin2//nfreq

    spectrum_var = res.get_spectrum_varname(point)

    # Getting record from time if given
    if time is not None:
        record = res.get_closest_record(time)

    # Reshaping to match nfreq*ntheta
    tmp_data = res.get_data_value(spectrum_var, record)\
        .reshape((nfreq, ntheta))

    # Adding frequencies as first column
    data = np.column_stack((freqs, tmp_data))

    res.close()
    # Building headers
    header = ['theta']

    # Defining if we are in radian or degree
    if radian:
        val = 2*np.pi
    else:
        val = 360.

    # Building angles array
    for i in range(ntheta):
        header.append(str(i*val/ntheta))

    return header, data
