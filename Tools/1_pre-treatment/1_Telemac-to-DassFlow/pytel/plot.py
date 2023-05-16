#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   @brief Simple plots
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from argparse import ArgumentParser
import numpy as np
from postel.plot_actions import \
        add_options_var, plot_var,\
        add_options_mesh2d, plot_mesh2d,\
        add_options_vertical_slice, plot_vertical_slice, \
        add_options_horizontal_slice, plot_horizontal_slice,\
        add_options_3d_scalar_map, \
        plot_timeseries_on_polyline, add_options_timeseries_on_polyline, \
        add_options_history, plot_history, \
        add_options_spe_ang, plot_spe_ang, \
        add_options_spe_freq, plot_spe_freq, \
        add_options_spe, plot_spe, \
        add_options_poly, plot_poly
from postel.plot_vnv import vnv_plot3d
from data_manip.extraction.shapefile_reader import read_shape_data
from data_manip.extraction.telemac_file import TelemacFile
from vvytel.report_class import Report


def add_options_report(subparser):
    """
    Defines options for var action

    @param subparser (ArgumentParser) The subparser

    @returns the update subparser
    """
    parser = subparser.add_parser('report',\
            help='Plot data from validation report')
    parser.add_argument("report_file", default=None, \
        help="Name of the report file")

    return subparser

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

def main():
    """
    Main function of the converter
    """
# ~~~~ Defines arguments
    parser = ArgumentParser()
    subparser = parser.add_subparsers(help='plot command to do', dest='command')

    subparser = add_options_mesh2d(subparser)
    subparser = add_options_var(subparser)
    subparser = add_options_3d_scalar_map(subparser)
    subparser = add_options_vertical_slice(subparser)
    subparser = add_options_horizontal_slice(subparser)
    subparser = add_options_timeseries_on_polyline(subparser)
    subparser = add_options_history(subparser)
    subparser = add_options_spe(subparser)
    subparser = add_options_spe_freq(subparser)
    subparser = add_options_spe_ang(subparser)
    subparser = add_options_report(subparser)
    subparser = add_options_poly(subparser)


    options = parser.parse_args()

    command = options.command

#   Actions to run
    if command == 'mesh2d':
        res = TelemacFile(options.input_file, bnd_file=options.bnd_file)
        plot_mesh2d(res, options.bnd,
                    options.liq_bnd, options.fig_name)
        res.close()
    elif command == 'var':
        res = TelemacFile(options.input_file)
        plot_var(res, options.var, options.record,
                 options.time, options.mesh, options.fig_name)
        res.close()
    elif command == '3d_scalar_map':
        res = TelemacFile(options.input_file)
        vnv_plot3d(\
                options.var, res, record=options.record,
                time=options.time, fig_name=options.fig_name)
        res.close()
    elif command == 'v_slice':
        res = TelemacFile(options.input_file)
        plot_vertical_slice(\
                res, options.var,
                options.poly, record=options.record,
                time=options.time, add_mesh=options.mesh,
                fig_name=options.fig_name)
        res.close()
    elif command == 'h_slice':
        res = TelemacFile(options.input_file)
        plot_horizontal_slice(\
                res, options.var,
                options.plane, options.record,
                options.time, options.mesh, options.fig_name)
        res.close()
    elif command == 'time_poly':
        res = TelemacFile(options.input_file)
        plot_timeseries_on_polyline(\
                res, options.var, options.poly,
                fig_name=options.fig_name)
        res.close()
    elif command == 'history':
        res = TelemacFile(options.input_file)
        plot_history(res, options.var, options.points,
                     fig_name=options.fig_name)
        res.close()
    elif command == 'spec':
        res = TelemacFile(options.input_file)
        plot_spe(res, options.point, options.record, options.time,
                 options.fig_name)
        res.close()

    elif command == 'spec-freq':
        res = TelemacFile(options.input_file)
        plot_spe_freq(res, options.points, options.record, options.time,
                      options.fig_name)
        res.close()
    elif command == 'spec-ang':
        res = TelemacFile(options.input_file)
        plot_spe_ang(res, options.points, options.record, options.time,
                     options.fig_name)
        res.close()
    elif command == 'report':

        rep = Report('', 'examples')

        rep.read(file_name=options.report_file)

        rep.plot_stats()

        del rep

    elif command == 'poly':
        # TODO: Handle names ?
        polys = read_shape_data(options.shp_file)
        plot_poly(polys[0])
    else:
        parser.print_help()


    sys.exit(0)

if __name__ == "__main__":
    main()
