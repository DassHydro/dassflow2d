#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   @brief Run a converions of mesh files using stbtel
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from argparse import ArgumentParser
from data_manip.conversion.convert_ecmwf import ecmwf2srf, ecmwf2srf_parser
from data_manip.conversion.convert_gebco import gebco2srf, gebco2srf_parser
from data_manip.conversion.convert_hycom import hycom2srf, hycom2srf_parser
from data_manip.conversion.convert_kenue import kenue2shp, kenue2shp_parser
from data_manip.conversion.dat2vtu import convert_drogues_file_to_vtu, \
                                          dat2vtu_parser
from data_manip.conversion.stbtel_converter import run_converter, \
                                                   stbtel_converter_parser
from data_manip.conversion.sis2gaia import sis2gaia_parser, sis2gaia
from data_manip.extraction.shapefile_reader import write_poly_txt, \
                 write_poly_i2s, read_shape_data, write_shape, read_poly_txt
from data_manip.extraction.extract_ptravers_res_to_geoc import \
        extract_ptravers_res_to_geoc, extrac_ptraver_parser
from pretel.convert_listing_courlis import\
        convert_listing_courlis, convert_courlis_parser
from pretel.generate_atm import generate_atm, generate_atm_parser
from pretel.convert_to_bnd import generate_bnd, generate_bnd_parser
from pretel.stbtel_refine import run_refine, stbtel_refine_parser
from vvytel.xml2py import xml2py, xml2py_parser
from vvytel.report_class import Report

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Yoann Audouin"
__date__ = "$21-Sep-2012 16:51:09$"

def report2xls_parser(subparser):
    """
    Add arguments for a conversion of report2xls

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser('report2xls',
                                  help='Converting validation report csv'\
                                       'into Xls file')
    parser.add_argument(
        "report_file",
        help="path of the report file.")
    parser.add_argument(
        "xls_file",
        help="path of the xls file.")
    parser.add_argument(
        "--mode",
        dest="mode", default='create',
        choices=['create', 'append', 'insert'],
        help="mode of xrite (create for a new file, "
             "append to add a new column at the end, "
             "insert to add a new column at the beginning")
    parser.add_argument(
        "--title",
        dest="title", default=None,
        help='Title give to the worksheet')
    parser.add_argument(
        "--job-id",
        dest="job_id", default=None,
        help='Name of the column of data')
    parser.add_argument(
        "--date",
        dest="date", default=None,
        help='Date of the column of data')
    parser.add_argument(
        "--guess",
        dest="guess", action='store_true', default=False,
        help='Will guess job_id, title and date from report_file name')
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose", action='store_true', default=False,
        help='Will print more information')
    parser.add_argument(
        "--max-report",
        dest="max_report", default=0, type=int,
        help='Fix max of report within the file (delete oldest one)')

    return subparser

def shp2i2s_parser(subparser):
    """
    Add arguments for a conversion of shp2i2s

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser(\
        'shp2i2s',
        help='Converting Shape file into i2s (BlueKenue polygon format)')
    parser.add_argument(
        "shape_file",
        help="path of the shape file.")
    parser.add_argument(
        "i2s_file",
        help="path of the i2s file.")
    parser.add_argument(
        "--poly-val", default=None, type=int, nargs="+",
        help="value for each polygon If only one is given it will be apllied "
             " to all the polygon.")

    return subparser

def shp2txt_parser(subparser):
    """
    Add arguments for a conversion of shp2txt

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser(\
        'shp2txt',
        help='Converting Shape file into txt (Telemac in_poly format)')
    parser.add_argument(
        "shape_file",
        help="path of the shape file.")
    parser.add_argument(
        "txt_file",
        help="path of the txt file.")

    return subparser

def txt2shp_parser(subparser):
    """
    Add arguments for a conversion of txt2shp

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser(\
        'txt2shp',
        help='Converting txt (Telemac in_poly format) into shape file')
    parser.add_argument(
        "txt_file",
        help="path of the txt file.")
    parser.add_argument(
        "shape_file",
        help="path of the shape file.")

    return subparser

def main():
    """
    Main function of the converter
    """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    parser = ArgumentParser()
    subparser = parser.add_subparsers(\
            help='converter command to do', dest='command')

    subparser = ecmwf2srf_parser(subparser)
    subparser = gebco2srf_parser(subparser)
    subparser = hycom2srf_parser(subparser)
    subparser = kenue2shp_parser(subparser)
    subparser = dat2vtu_parser(subparser)
    subparser = generate_atm_parser(subparser)
    subparser = generate_bnd_parser(subparser)
    subparser = stbtel_converter_parser(subparser, 'srf2med',\
            'Conversion from serafin (single or double precision) to MED')
    subparser = stbtel_converter_parser(subparser, 'srf2vtk',\
            'Conversion from serafin (single or double precision) to ascii VTK')
    subparser = stbtel_converter_parser(subparser, 'med2srf',\
            'Conversion from MED to serafin single precision')
    subparser = stbtel_converter_parser(subparser, 'med2srfd',\
            'Conversion from MED to serafin double precision')
    subparser = stbtel_refine_parser(subparser)
    subparser = sis2gaia_parser(subparser)
    subparser = xml2py_parser(subparser)
    subparser = extrac_ptraver_parser(subparser)
    subparser = convert_courlis_parser(subparser)
    subparser = report2xls_parser(subparser)
    subparser = shp2i2s_parser(subparser)
    subparser = shp2txt_parser(subparser)
    subparser = txt2shp_parser(subparser)

    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads code name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if options.command == 'ecmwf2srf':
        ecmwf2srf(options.tfrom, options.tstop, options.blcorner,
                  options.trcorner, options.dataset, options.stream,
                  options.root_name)
    elif options.command == 'gebco2srf':
        gebco2srf(options.gebco_file, options.abval, options.beval,
                  options.axp, options.sph2ll, options.ll2sph,
                  options.ll2utm, options.utm2ll)
    elif options.command == 'hycom2srf':
        hycom2srf(options.tfrom, options.tstop, options.blcorner,
                  options.trcorner, options.root_name, options.t2d)
    elif options.command == 'kenue2shp':
        kenue2shp(options.args)
    elif options.command == 'dat2vtu':
        convert_drogues_file_to_vtu(options.drogues_file, options.vtu_file)
    elif options.command == 'generate_atm':
        generate_atm(options.geo_file, options.slf_file, options.atm_file,
                     options.ll2utm)
    elif options.command == 'generate_bnd':
        varnames = []
        for var in options.varnames.split(";"):
            if len(var) < 16:
                name = var + ' '*(16-len(var))
            else:
                name = var[:16]
            varnames.append(name)

        varunits = []
        for var in options.varunits.split(";"):
            if len(var) < 16:
                name = var + ' '*(16-len(var))
            else:
                name = var[:16]
            varunits.append(name)

        generate_bnd(options.cli_file, options.geo_file, options.slf_file,
                     options.bnd_file,
                     varnames,
                     varunits)
    elif options.command in ['srf2med', 'srf2vtk', 'med2srf', 'med2srfd']:
        run_converter(options.command, options.input_file, options.output_file,
                      options.root_dir, options.bnd_file, options.log_file,
                      options.ndomains, options.srf_bnd, options.debug)
    elif options.command == "refine":
        run_refine(options.input_file, options.output_file, options.root_dir,
                   options.bnd_file)
    elif options.command == "sis2gaia":
        sis2gaia(options.sis_cas, options.gaia_cas)
    elif options.command == "xml2py":
        xml2py(options.input_file, options.out_file, options.skip)
    elif options.command == "extract_ptravers_res_to_geoc":
        extract_ptravers_res_to_geoc(options.args, options.outputFileName,
                                     options.record)
    elif options.command == "convert_listing_courlis":
        convert_listing_courlis(options.args, options.outputFileName,
                                options.write_opt, options.budget,
                                options.txt_format, options.csv_format,
                                options.xlsx_format, options.time_check)
    elif options.command == "report2xls":
        rep = Report('', 'examples')

        if options.guess:
            import re
            test = re.compile(r"(?P<name>[-a-zA-Z]+)_(?P<job_id>[0-9]+)_(?P<config>[-.a-zA-Z0-9]+)_[-.a-zA-Z]+_[-.a-zA-Z]+_(?P<date>[-0-9]+)-[0-9hmins]+.csv")
            proc = test.search(options.report_file)

            options.job_id = proc.group('job_id')
            options.date = proc.group('date')
            options.title = proc.group('config')

        rep.read(file_name=options.report_file)
        rep.write2xls(options.xls_file, options.mode, options.title,
                      options.job_id, options.date, verbose=options.verbose,
                      max_report=options.max_report)
    elif options.command == "shp2txt":
        poly = read_shape_data(options.shape_file)
        write_poly_txt(poly, options.txt_file)
        print(" ~> Created {}".format(options.txt_file))
    elif options.command == "shp2i2s":
        poly = read_shape_data(options.shape_file)
        poly_val = None
        if options.poly_val is not None:
            print(options.poly_val)
            if len(options.poly_val) == 1:
                poly_val = [options.poly_val[0]]*len(poly)
            else:
                poly_val = options.poly_val
        write_poly_i2s(poly, options.i2s_file, poly_val=poly_val)
        print(" ~> Created {}".format(options.i2s_file))
    elif options.command == 'txt2shp':
        polys = read_poly_txt(options.txt_file)
        write_shape(options.shape_file, polys)
    else:
        parser.print_help()

    sys.exit(0)

if __name__ == "__main__":
    main()
