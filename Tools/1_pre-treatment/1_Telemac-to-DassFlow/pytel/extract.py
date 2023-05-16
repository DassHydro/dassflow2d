#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   @brief Simple extraction from csv files
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import path
from argparse import ArgumentParser
import numpy as np
from data_manip.extraction.extract_actions import \
       extract_timeseries, add_options_timeseries, \
       extract_mesh2d, add_options_mesh2d, \
       extract_spectrum, add_options_spectrum


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

    subparser = add_options_timeseries(subparser)
    subparser = add_options_mesh2d(subparser)
    subparser = add_options_spectrum(subparser)

    options = parser.parse_args()

    command = options.command

#   Actions to run
    if command == 'timeseries':
        variables = options.var.split(',')
        for var in variables:
            header, data = extract_timeseries(options.file_name, var,
                                              nodes=options.nodes,
                                              points=options.points)
            # Writting csv file
            csv_name = options.csv_name
            if len(variables) > 1:
                root, ext = path.splitext(options.csv_name)
                csv_name = "{}_{}{}".format(root, var.strip(' ').replace(' ', '_'), ext)

            np.savetxt(csv_name, data, header=options.delimiter.join(header),
                       delimiter=options.delimiter)

    elif command == "mesh2d":
        header, data = extract_mesh2d(options.file_name)
        # Writting csv file
        np.savetxt(options.csv_name, data, header=options.delimiter.join(header),
                   delimiter=options.delimiter)

    elif command == 'spectrum':
        header, data = extract_spectrum(\
                options.file_name, options.point,
                radian=options.radian,
                time=options.time, record=options.record)
        # Writting csv file
        np.savetxt(options.csv_name, data, header=options.delimiter.join(header),
                   delimiter=options.delimiter)
    else:
        parser.print_help()


    sys.exit(0)

if __name__ == "__main__":
    main()
