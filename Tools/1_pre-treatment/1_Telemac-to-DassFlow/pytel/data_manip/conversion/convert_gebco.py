r"""@author TELEMAC-MASCARET Consortium

     @brief
     Tools for handling conversions to-from GEBCO server files

     @details
     Contains server read functions to convert to SELAFIN file
"""

from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from os import path
import numpy as np
from matplotlib.tri import Triangulation
# ~~> dependencies towards other modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
import data_manip.conversion.convert_utm as utm
from pretel import polygons

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#


class Gebco(Selafin):

    def __init__(self, fname, vals=(None, None)):

        """
        @brief Main class that fetches the GEBCO bathymetry data and puts it
                 into the SELAFIN data format

        @param fname (string) the name of the GEBCO file
        @param vals (float, float) z min, max
        """

        # ~~> empty SELAFIN
        Selafin.__init__(self, '')
        self.datetime = []

        # ~~> variables
        self.title = ''
        self.nbv1 = 1  # bathymetry only
        self.nvar = self.nbv1
        self.varindex = range(self.nvar)
        self.varnames = ['BOTTOM          ']
        self.varunits = ['M               ']

        print('     +> header')
        # ~~> load header (ASC type)
        gebcofile = open(fname, 'r')
        # ~~
        gline = []
        gline.append(gebcofile.readline().split())
        if gline[-1][0].lower() == "ncols":
            nx1d = int(gline[-1][1])
        else:
            raise TelemacException(
                    '.. Could not read this file format. '
                    'Key ncols expected here.')
        gline.append(gebcofile.readline().split())
        if gline[-1][0].lower() == "nrows":
            ny1d = int(gline[-1][1])
        else:
            raise TelemacException(
                    '.. Could not read this file format. '
                    'Key nrows expected here.')
        gline.append(gebcofile.readline().split())
        if gline[-1][0].lower() == "xllcorner":
            xllcorner = np.float(gline[-1][1])
        else:
            raise TelemacException(
                    '.. Could not read this file format. '
                    'Key xllcorner expected here.')
        gline.append(gebcofile.readline().split())
        if gline[-1][0].lower() == "yllcorner":
            yllcorner = np.float(gline[-1][1])
        else:
            raise TelemacException(
                    '.. Could not read this file format. '
                    'Key yllcorner expected here.')
        gline.append(gebcofile.readline().split())
        if gline[-1][0].lower() == "cellsize":
            xdim = np.float(gline[-1][1])
            ydim = xdim
        elif gline[-1][0].lower() in ["xdim", "dx"]:
            xdim = np.float(gline[-1][1])
            gline.append(gebcofile.readline().split())
            if gline[-1][0].lower() in ["ydim", "dy"]:
                ydim = np.float(gline[-1][1])
            else:
                raise TelemacException(
                        '.. Could not read this file format.'
                        ' Key ydim expected here.')
        else:
            raise TelemacException(
                    '.. Could not read this file format. '
                    'Key cellsize or xdim expected here.')
        gline.append(gebcofile.readline().split())
        if gline[-1][0].lower() == "nodata_value":
            nodata_value = int(gline[-1][1])
        else:
            raise TelemacException(
                    '.. Could not read this file format. '
                    'Key NODATA_value expected here.')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        gebcofile.close()

        print('     +> bathymetry')
        # ~~> load ASCII content, ignoring the header lines
        z = np.loadtxt(fname, skiprows=len(gline)).T.ravel()
        print('     +> filtered connectivity')
        # ~~> temporary ikle
        aval = min(z) - 1
        if vals[0] is not None:
            aval = float(vals[0])
        bval = max(z) + 1
        if vals[1] is not None:
            bval = float(vals[1])
        ielem = 0
        pbar = ProgressBar(maxval=2*(nx1d-1)*(ny1d-1)).start()
        ikle3 = - np.ones((2*(nx1d-1)*(ny1d-1), 3), dtype=np.int)
        for i in range(1, nx1d):
            for j in range(1, ny1d):
                ipoin = (i-1)*ny1d + j - 1
                # ~~> first triangle
                if (aval < z[ipoin] < bval) and \
                   (aval < z[ipoin + ny1d] < bval) and \
                   (aval < z[ipoin + 1] < bval):
                    ikle3[ielem] = [ipoin, ipoin + 1, ipoin + ny1d]
                ielem = ielem + 1
                pbar.update(ielem)
                # ~~> second triangle
                if (aval < z[ipoin + ny1d] < bval) and \
                   (aval < z[ipoin + ny1d + 1] < bval) and \
                   (aval < z[ipoin + 1] < bval):
                    ikle3[ielem] = [ipoin + ny1d, ipoin + 1, ipoin + ny1d + 1]
                ielem = ielem + 1
                pbar.update(ielem)
        pbar.finish()

        print('     +> renumbered connectivity')
        # ~~> intermediate connectivity
        gikle = ikle3[np.not_equal(*(np.sort(ikle3).T[0::2]))]
        knolg = np.unique(np.ravel(gikle))
        knogl = dict(zip(knolg, range(len(knolg))))
        # ~~> final connectivity
        self.ikle3 = - np.ones_like(gikle, dtype=np.int)
        pbar = ProgressBar(maxval=len(gikle)).start()
        for k in range(len(gikle)):
            self.ikle3[k] = [knogl[gikle[k][0]],
                             knogl[gikle[k][1]],
                             knogl[gikle[k][2]]]
            pbar.update(k)
        pbar.finish()

        print('     +> mesh x,y,z')
        # ~~> defines grid
        x = xllcorner + xdim * np.arange(nx1d, dtype=np.float) - xdim / 2.
        y = yllcorner - ydim * np.arange(ny1d, dtype=np.float) + \
            ydim * ny1d - ydim / 2.
        self.meshx = np.tile(x, ny1d).reshape(ny1d, nx1d).T.ravel()[knolg]
        self.meshy = np.tile(y, nx1d)[knolg]
        self.z = z[knolg]

        print('     +> sizes')
        # ~~> sizes
        self.nplan = 1
        self.ndp2 = 3
        self.ndp3 = self.ndp2
        self.npoin2 = len(self.meshx)
        self.npoin3 = self.npoin2
        self.nelem2 = len(self.ikle3)
        self.nelem3 = self.nelem2
        self.iparam = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

        print('     +> boundaries')
        # ~~> establish neighborhood
        neighbours = Triangulation(self.meshx, self.meshy, self.ikle3)\
            .get_cpp_triangulation().get_neighbors()
        # ~~> build the enssemble of boundary segments
        ebounds = []
        print('        - identify')
        pbar = ProgressBar(maxval=self.nelem3).start()
        for i in range(self.nelem3):
            if neighbours[i, 0] < 0:
                ebounds.append([self.ikle3[i][0], self.ikle3[i][1]])
            if neighbours[i, 1] < 0:
                ebounds.append([self.ikle3[i][1], self.ikle3[i][2]])
            if neighbours[i, 2] < 0:
                ebounds.append([self.ikle3[i][2], self.ikle3[i][0]])
            pbar.update(i)
        pbar.finish()
        # ~~> assemble the enssemble of boundary segments
        print('        - assemble')
        pbounds = polygons.join_segments(ebounds)
        # ~~> define ipobO from an arbitrary start point
        print('        - set')
        self.ipob3 = np.zeros(self.npoin3, dtype=np.int)
        iptfr = 0
        for bound in pbounds:
            for n in bound[1:]:
                iptfr += 1
                self.ipob3[n] = iptfr
        self.ipob2 = self.ipob3

    def put_content(self, file_name, showbar=True):

        """
        @brief Writes the GEBCO file to the SELAFIN format

        @param file_name (string): name of the SELAFIN file to write
        @param showbar (boolean): weather to show progress bar
        """

        # ~~> new SELAFIN writer
        self.fole = {}
        self.fole.update({'hook': open(file_name, 'wb')})
        self.fole.update({'name': file_name})
        self.fole.update({'endian': ">"})     # big endian
        self.fole.update({'float': ('f', 4)})  # single precision

        print('     +> Write SELAFIN header')
        self.append_header_slf()

        print('     +> Write SELAFIN core')
        self.append_core_time_slf(0.0)
        self.append_core_vars_slf([self.z])
        self.fole['hook'].close()


def gebco2srf(gebco_file, abval, beval, axp, sph2ll, ll2sph, ll2utm, utm2ll):
    """
    """
    if not path.exists(gebco_file):
        raise TelemacException(
                '... the provided GEBCO file does '
                'not seem to exist: {}\n\n'.format(gebco_file))
    head, _ = path.splitext(gebco_file)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Download the GEBCO file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\n'+72*'~'+'\n')
    print('\nLoading the GEBCO file\n')
    gebco2slf = Gebco(gebco_file, (abval, beval))

    gebco2slf.meshx = gebco2slf.meshx + float(axp)
    if sph2ll is not None:
        radius = 6371000.
        long0, lat0 = sph2ll.split(":")
        long0 = np.deg2rad(float(long0))
        lat0 = np.deg2rad(float(lat0))
        const = np.tan(lat0/2. + np.pi/4.)
        gebco2slf.meshx = np.rad2deg(gebco2slf.meshx/radius + long0)
        tmp = np.arctan(const*np.exp(gebco2slf.meshy/radius))
        gebco2slf.meshy = np.rad2deg(2.*tmp - np.pi/2.)
    if ll2sph is not None:
        radius = 6371000.
        long0, lat0 = ll2sph.split(":")
        long0 = np.deg2rad(float(long0))
        lat0 = np.deg2rad(float(lat0))
        gebco2slf.meshx = radius * (np.deg2rad(gebco2slf.meshx) - long0)
        gebco2slf.meshy = \
            radius * (np.log(np.tan(np.deg2rad(gebco2slf.meshy)/2. + np.pi/4.))
                      - np.log(np.tan(lat0/2. + np.pi/4.)))
    if ll2utm is not None:
        zone = int(ll2utm)
        gebco2slf.meshx, gebco2slf.meshy, zone = \
            utm.from_lat_long(gebco2slf.meshx, gebco2slf.meshy, zone)
    if utm2ll is not None:
        zone = int(utm2ll)
        gebco2slf.meshx, gebco2slf.meshy = \
            utm.to_lat_long(gebco2slf.meshx, gebco2slf.meshy, zone)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Convert to SELAFIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\n'+72*'~'+'\n')
    print('\nConverting into SELAFIN\n')
    gebco2slf.put_content(head+'.slf')


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
def gebco2srf_parser(subparser):
    """
    Add arguments for a gebco2srf conversion

    @param subparser (argumentParser) argument parser

    @return (argumentParser) the updated argument parser
    """
    parser = subparser.add_parser('gebco2srf',
                                  help='Extracting GEBCO data into\
                                        a serafin-format file')
    parser.add_argument(
        "gebco_file", default='',
        help="Name of the gebco file to extract data from.")
    parser.add_argument(
        "--above",
        dest="abval", default=None,
        help="select only the values above")
    parser.add_argument(
        "--below",
        dest="beval", default=None,
        help="select only the values below")
    parser.add_argument(
        "--sph2ll",
        dest="sph2ll", default=None,
        help="convert from spherical to longitude-latitude")
    parser.add_argument(
        "--ll2sph",
        dest="ll2sph", default=None,
        help="convert from longitude-latitude to spherical")
    parser.add_argument(
        "--ll2utm",
        dest="ll2utm", default=None,
        help="convert from longitude-latitude to UTM")
    parser.add_argument(
        "--utm2ll",
        dest="utm2ll", default=None,
        help="convert from UTM to longitude-latitude")
    parser.add_argument(
        "--X+?",
        dest="axp", default="0",
        help="adds to the meshx")

    return subparser


__author__ = "Juliette Parisi"
__date__ = "$12-Oct-2015 08:51:29$"
