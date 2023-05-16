#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   Manipulation of Telemac files (mesh, results files)
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from argparse import ArgumentParser
# ~~> dependencies towards other modules
from pretel.compute_weight import connect_tel2tom
from pretel.extract_contour import extract_contour, write_gis_file
from pretel.clc import add_clc_in_file
from pretel.manip_telfile import scan, diff
from utils.exceptions import TelemacException

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

def tel2tom_parser(subparser):
    """ Generate parser for tel2tom """
    parser = subparser.add_parser('tel2tom',\
            help="Compute weight for each node and index to interpolate from "
                 "telemac2d file to tomawac file and reverse")
    parser.add_argument("t2d_file", metavar="Telemac File", help="First file")
    parser.add_argument("tom_file", metavar="Tomawac File", help="Second file")
    parser.add_argument("--t2d-contour", default=None,
                        help="Polygon to apply to Telemac2d file")
    parser.add_argument("--tom-contour", default=None,
                        help="Polygon to apply to Tomawac file")

    return subparser

def contour_parser(subparser):
    """ Generate parser for extract_contour """
    parser = subparser.add_parser('contour',\
            help="Compute contour (in Shape file format)from mesh file")
    parser.add_argument("mesh_file", metavar="Telemac File", help="Mesh file")
    parser.add_argument("shp_file", help="Shape file")
    parser.add_argument("-b", "--bnd-file", default=None, help="Boundary file")

    return subparser

def clc_parser(subparser):
    """ Generate parser for extract_contour """
    parser = subparser.add_parser('clc',\
            help="Compute Corinne Land Cover data on mesh")
    parser.add_argument("mesh_file", metavar="Telemac File", help="Mesh file")
    parser.add_argument(
        "-v", dest="varname", default="FRICTION",
        help="Name of the variable containing Corinne Land Cover data "
             "(default: FRICTION)")
    parser.add_argument("--year", default=2012,
                        choices=[2012],
                        help="Year to extract data from")
    parser.add_argument("-o", dest="res_file", metavar="Telemac File",
                        help="Output file if not given will write in mesh_file")

    return subparser

def scan_parser(subparser):
    """ Generate parser for scan """
    parser = subparser.add_parser('scan',\
            help="Dumping file information")
    parser.add_argument("tel_file", metavar="Telemac File", help="File to scan")
    parser.add_argument(
        "-b", "--boundary_file", dest="bnd_file",
        metavar="Telemac Boundary File", help="Boundary file")
    parser.add_argument(
        "--data", default=False, action='store_true',
        help="Display data information")

    return subparser

def diff_parser(subparser):
    """ Generate parser for diff """
    parser = subparser.add_parser('diff',\
        help="Creates a file containing the difference betwenn two other files")
    parser.add_argument("tel_file1", help="File1")
    parser.add_argument("tel_file2", help="File2")
    parser.add_argument("diff_file", help="File2-file1")

    return subparser


def main():
    """
    Main function
    """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+'~'*72+'\n')

    parser = ArgumentParser()
    subparser = parser.add_subparsers(\
            help='run_selafin commands to do', dest='command')

    subparser = tel2tom_parser(subparser)
    subparser = contour_parser(subparser)
    subparser = clc_parser(subparser)
    subparser = scan_parser(subparser)
    subparser = diff_parser(subparser)

    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads code name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if options.command == 'tel2tom':
        connect_tel2tom(options.t2d_file, options.tom_file,
                        contour_tom=options.tom_contour,
                        contour_tel=options.t2d_contour)
    elif options.command == 'contour':
        domains_bnd = extract_contour(options.mesh_file, \
                        bnd_file=options.bnd_file)

        write_gis_file(domains_bnd, options.shp_file)
    elif options.command == 'clc':
        add_clc_in_file(options.mesh_file, options.varname, options.year,
                        options.res_file)
    elif options.command == 'scan':
        scan(options.tel_file, options.bnd_file, options.data)
    elif options.command == 'diff':
        diff(options.tel_file1, options.tel_file2, options.diff_file)
    else:
        raise TelemacException(\
                '\nDo not know what to do with '
                'this code name: {}'.format(options.command))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)

if __name__ == "__main__":
    main()
