r"""@author TELEMAC-MASCARET Consortium

    @brief
    Tools for handling conversions to-from ECMWF server files

    @details
    Contains server read functions to convert to SELAFIN file
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from os import path, listdir
from ast import literal_eval
# ~~> dependencies towards other pytel scripts
# ~~> dependencies towards other modules
from data_manip.formats.grib import Grib
from data_manip.formats.ecmwf import Ecmwf
from utils.exceptions import TelemacException

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

try:
    import pygrib
    IMPORTGRIB = True
except ImportError:
    IMPORTGRIB = False

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban"
__date__ = "$12-Dec-2014 08:51:29$"


def build_period(tfrom, tstop):
    """
    Define the time period for the data request

    @param tfrom (string) First date (ex "1972-13-07")
    @param tstop (string) Last date (ex "1980-12-31")
    """
    period = [[], []]
    if tfrom is not None:
        for i in tfrom.split('-'):
            period[0].append(i)
    else:
        raise TelemacException(
                '... could not find your from date. '
                'Please use --from option '
                '(- delimited, no spaces).\n\n')
    if tstop is not None:
        for i in tstop.split('-'):
            period[1].append(i)
    else:
        raise TelemacException(
                '... could not find your stop date. '
                'Please use --stop option '
                '(- delimited, no spaces).\n\n')
    return period


def build_box(blcorner, trcorner):
    """
    Define the box for the data request

    @param blcorner (string) Bottom Left cnorner
    @param trcorner (string) Top Right corner

    @return The model box
    """
    modelbox = [[], []]
    if blcorner is not None:
        x, y = literal_eval(blcorner)
        modelbox[0].append(str(x))
        modelbox[0].append(str(y))
    else:
        raise TelemacException(
                '... could not find your bounding box bottom left'
                'corner. Please use --bl option (, delimited, no '
                'spaces).\n\n')
    if trcorner is not None:
        x, y = literal_eval(trcorner)
        modelbox[1].append(str(x))
        modelbox[1].append(str(y))
    else:
        raise TelemacException(
                '... could not find your bounding box top right '
                'corner. Please use --tr option (, delimited, no '
                'spaces).\n\n')
    return modelbox


def build_data_request(modelbox, period, dataset, stream, res_name):
    """
    Build the dictionnary for the data request

    @param modelbox (List) Area to extract
    @param period (List) Time period to extract
    @param dataset (string) Type of dataset to use
    @param stream (string) Type of stream
    @param res_name (string) Name of the result file

    @return The dictionary containing the request information and the filetype
    """
    head, _ = path.splitext(res_name)

    if dataset == 'interim':
        download = True
        req = {
            'dataset': "interim",
            'step': "0",
            'number': "all",
            'levtype': "sfc",      # surface or single level
            'date': "2011-02-15/to/2011-06-15",
            'time': "00/06/12/18",
            # 'origin'  : "all",
            'type': "an",       # analysis
            # because these are ambigous: "ap/u10/v10/t2m"
            'param': "134.128/165.128/166.128/167.128",
            'area': "41/140/34/147",
            'grid': "0.125/0.125",
            'target': head+'.nc',
            'format': 'netcdf',
            'class': "ei",
            'stream': stream
              }
        if stream == 'wave':
            req['param'] = "swh/mwp/mwd"
        # only supported file type
        file_type = 'netcdf'
    elif dataset == 'era5':
        download = True
        req = {
            'dataset': "era5",
            'step': "0",
            'number': "all",
            'levtype': "sfc",      # surface or single level
            'date': "2011-02-15/to/2011-06-15",
            'time': "00/01/02/03/04/05/06/07/08/09/10/11/12/13/"
                    "14/15/16/17/18/19/20/21/22/23",
            # 'origin'  : "all",
            'type': "an",       # analysis
            # because these are ambigous: "ap/u10/v10/t2m"
            'param': "134.128/165.128/166.128/167.128",
            'area': "41/140/34/147",
            'grid': "0.125/0.125",
            'target': head+'.nc',
            'format': 'netcdf',
            'class': "ea",
            'stream': stream
              }
        if stream == 'wave':
            req['param'] = "swh/mwp/mwd"
        if stream == 'spec':
            req['param'] = "2dfd"
            req['stream'] = 'wave'
        # only supported file type
        file_type = 'netcdf'
    else:
        download = False
        # ~~> GRIB files within a directory
        if path.isdir(dataset):
            # only supported file type
            file_type = 'grib'
            dataset = [path.join(dataset, d)
                       for d in listdir(dataset)]
        # ~~> actual file
        elif path.exists(dataset):
            if path.splitext(dataset)[1] == '.nc':
                file_type = 'netcdf'
            else:
                file_type = 'grib'
            dataset = [dataset]
        else:
            raise TelemacException(
                    '... could not find your local files.\n\n')
        # ~~> In case of GRIB files
        if file_type == 'grib':
            if not IMPORTGRIB:
                raise TelemacException(
                    '... could not import pygrib on this platform.\n\n')
        # ~~> In case of GRIB files
        req = {
            'date': "2011-02-15/to/2011-06-15",
            'time': "00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15"
                    "/16/17/18/19/20/21/22/23",
            # because these are ambigous: "ap/u10/v10"
            'param': "151/165/166",
            'area': "41/140/34/147",
            'grid': "0.125/0.125",
            'target': dataset[0]
              }
        if stream == 'spec':
            req['param'] = "251"

    req['date'] = '-'.join(period[0]) + '/to/' + '-'.join(period[1])
    req['area'] = '/'.join(modelbox[0]) + '/' + \
                  '/'.join(modelbox[1])
    return req, file_type, download


def write_file(req, file_type, download, dataset, stream, period, root_name):
    """
    Writing the Serafin file according to data request

    @param req (dict) Contains all the request data
    @param file_type (string) Format on input data
    @param download (boolean) wether or not we need to dwnload stuff
    @param dataset (string) Type of dataset
    @param stream (string) Type of stream
    @param period (List) Time period
    @param root_name (string) Name of the result file
    """
# ~~~~ Loading up the GRIB file~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
    head, _ = path.splitext(root_name)

    if file_type == 'grib':

        if download:
            raise TelemacException(
                    '... I am not programmed to '
                    'download grib files directly.\n\n')
        print('\n\n'+72*'~'+'\n')
        print('\nLoading essentials from the GRIB\n')
        grb2slf = Grib(dataset, req, stream)

        grb2slf.set_geometry()

        if stream == 'spec':
            print('\n\n'+72*'~'+'\n')
            print('\nSpecial case for spectral file\n')
            grb2slf.put_geometry('geo_'+head+'.slf')
            grb2slf.set_spectral()

        print('\n\n'+72*'~'+'\n')
        print('\nConverting grib file(s) into SELAFIN\n')
        grb2slf.put_content(root_name)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Downloading the NetCDF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Unfortunately, I did not manage to access the NetCDF file remotely
    elif file_type == 'netcdf':

        ecmwf2slf = Ecmwf(period, req)
        if download:
            print('\n\n'+72*'~'+'\n')
            print('\nMaking an ECMWF request\n')
            ecmwf2slf.connect_to_ecmwf("datasets/%s" % (req['dataset']))

            print('\n\n'+72*'~'+'\n')
            print('\nHaving to download the ECMWF file first\n')
            ecmwf2slf.download_ecmwf()
            print("   ~> download completed.")

        ecmwf2slf.open_ecmwf()
        ecmwf2slf.set_geometry()

        if stream == 'spec':
            print('\n\n'+72*'~'+'\n')
            print('\nSpecial case for spectral file\n')
            ecmwf2slf.put_geometry('geo_'+head+'.slf')
            ecmwf2slf.set_spectral()

        print('\n\n'+72*'~'+'\n')
        print('\nConverting netcdf file into SELAFIN\n')
        ecmwf2slf.put_content(root_name, stream)


def ecmwf2srf_parser(subparser):
    """
    Add arguments for a cemwf2srf conversion

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser('ecmwf2srf',
                                  help='extracting ECMWF data into\
                                        a serafin-format file')
    parser.add_argument(
        "root_name", default='',
        help="specify the root name of the resulting SELAFIN file.")
    parser.add_argument(
        "-f", "--from",
        dest="tfrom", default=None,
        help='specify the first date included. Example: "1972-13-07"')
    parser.add_argument(
        "-s", "--stop",
        dest="tstop", default=None,
        help='specify the last date included "1980-12-31"')
    parser.add_argument(
        "--bl", dest="blcorner", default=None,
        help='specify the bottom left corner. Example: "(25,-117)"')
    parser.add_argument(
        "--tr", dest="trcorner", default=None,
        help='specify the top right corner. Example: "(27,-110)"')
    parser.add_argument(
        "--dataset",
        dest="dataset", default='interim',
        help="type of dataset requested either 'interim', 'era5' or the name "
             "of a grib or netcdf file or dirctory contaiing grib or netcdf "
             "files, depending on periods and resolution, etc.")
    parser.add_argument(
        "--stream",
        dest="stream", default='oper',
        choices=set(("oper", "wave", "spec")),
        help="type of stream requested either 'oper' (atmospheric) or "
             "'wave' (waves), etc.")

    return subparser


def ecmwf2srf(tfrom, tstop, blcorner, trcorner, dataset, stream, root_name):
    """ Main function of conxvertECMWF """

    # Arbitrary 6-day period
    period = build_period(tfrom, tstop)

    # arbitrary box (small pieve of the atlantic side of Mexico)
    modelbox = build_box(blcorner, trcorner)

    req, file_type, download = build_data_request(modelbox, period,
                                                  dataset,
                                                  stream,
                                                  root_name)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    write_file(req, file_type, download, dataset, stream,
               period, root_name)
