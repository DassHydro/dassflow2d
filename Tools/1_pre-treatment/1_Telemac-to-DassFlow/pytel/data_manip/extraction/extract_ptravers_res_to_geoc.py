#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author TELEMAC-MASCARET Consortium
@brief Extract from ptravers Courlis res file to Courlis geometry file at a
given record
"""

import argparse
from data_manip.formats.mascaret_file import MascaretFile, Reach
from data_manip.formats.mascaretgeo_file import MascaretGeoFile


def extract_ptravers_res_to_geoc(name, output_name, record):
    """
    Extract from ptravers Courlis res file to Croulise geometry file at a given
    record

    @param name (str) Name of the file from which to extract
    @param output_name (str) Name of the output file
    @param record (int) Record to extract
    """

    courlis_ptravers = MascaretFile(name, fformat='ptravers')
    courlis_geo = MascaretGeoFile(output_name, mode='write')

    if output_name is '_prev.geoC':
        short_name = ".".join(name.split(".")[:-1])
        output_name = '{}_time_{}s.geoC'.format(short_name,
                                                courlis_ptravers.times[record])

    print('\nExtracting time {}s from {} to {}\n'
          .format(courlis_ptravers.times[record], name, output_name))

    # minus 4 because there are 4 variables in addition to layers elevation
    # minus 1 because it is elevation of layer, river bottom and hard bottom,
    # not thicknesses of layers
    nlayers = courlis_ptravers.nsectionvar - 4 - 1

    section_vars_indexes = [i for i in range(0, nlayers + 2)]

    _, layers_elevation = courlis_ptravers.get_values(
                             record,
                             get_section_values=True,
                             section_vars_indexes=section_vars_indexes)

    # Only one reach in Courlis
    reach = Reach(1, "Reach")
    courlis_geo.has_layers = True
    courlis_geo.nlayers = nlayers
    for i in range(nlayers+1):
        layer = courlis_ptravers.section_varnames_dict['names'][i+1]
        courlis_geo.layer_names.append(layer)

    for i in range(courlis_ptravers.nsections):

        section = courlis_ptravers.reaches[1].sections[i+1]
        dist = layers_elevation[1][i][0]
        z_list = layers_elevation[1][i][1]
        section.set_points_from_trans(dist, z_list)

        for j in range(nlayers):
            thickness_table = layers_elevation[1][i][j+1] - \
                              layers_elevation[1][i][j+2]
            section.add_layer(thickness_table)

        reach.add_section(section)

    courlis_geo.add_reach(reach)

    courlis_geo.save(output_name)


def extrac_ptraver_parser(subparser):
    """
    Build subparser

    @param subparser (ArgumentParser) The parser

    @returns (ArgumentParser) The updated parser
    """
    parser = subparser.add_parser('extract_ptravers_res_to_geoc',
                                  help='Extract from ptravers Courlis\
                                   res file to Courlis geometry'
                                  'file at a given record')
    parser.add_argument("args", metavar='ptravers Courlis file')
    parser.add_argument("-o", "--output",
                        dest="outputFileName",
                        default='_prev.geoC',
                        help="Option to give a name to the output Courlis "
                             "geometry file")
    parser.add_argument("-r", "--time-record",
                        dest="record", type=int, default='-1',
                        help="Record number of the state you want\
                              to extract from "
                             "ptravers to a geometry courlis file")

    return subparser


def main():
    """
    Main function
    """

    parser = argparse.ArgumentParser(
        description='Extract from ptravers Courlis res file to\
                        Courlis geometry'
                    'file at a given record')

    parser.add_argument("args", metavar='ptravers Courlis file')
    parser.add_argument("-o", "--output",
                        dest="outputFileName",
                        default='_prev.geoC',
                        help="Option to give a name to the output Courlis "
                             "geometry file")
    parser.add_argument("-r", "--time-record",
                        dest="record", type=int, default='-1',
                        help="Record number of the state you want\
                                to extract from "
                             "ptravers to a geometry courlis file")

    options = parser.parse_args()

    extract_ptravers_res_to_geoc(options.args, options.outputFileName,
                                 options.record)


if __name__ == '__main__':
    main()
