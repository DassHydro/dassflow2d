r"""@author Juliette C.E. Parisi and Sebastien E. Bourban

    @note ... this work is based on a collaborative effort between
  .________.                                                          ,--.
  |        |                                                      .  (  (
  |,-.    /   HR Wallingford                EDF - LNHE           / \_ \_/ .--.
  /   \  /    Howbery Park,                 6, quai Watier       \   )   /_   )
    ,.  `'     Wallingford, Oxfordshire      78401 Cedex           `-'_  __ `--
  /  \   /    OX10 8BA, United Kingdom      Chatou, France        __/ \ \ `.
 /    `-'|    www.hrwallingford.com         innovation.edf.com   |    )  )  )
!________!                                                        `--'   `--

    @brief
    Tools for handling conversions to-from HYCOM server files

    @details
    Contains server read functions to convert to SELAFIN file

    @history 30/04/2013 -- Sebastien E. Bourban and Juliette C.E. Parisi
    Complete re-write of the work prototype by Laure Grignon, including
    the production of 2D and 3D Selafin files.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> python 3 compatibility
# ~~> dependencies towards standard python
import sys
from os import path
import time
from datetime import datetime
import numpy as np
# ~~> dependencies towards other modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
# ~~> dependencies towards other pytel scripts
sys.path.append(path.join(path.dirname(sys.argv[0]), '..'))

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#
# Current Experiments
#    91.2 (18-Apr-2016 to Present)
#    => url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2'
#    91.1 (05-Apr-2014 to 18-Apr-2016)
#    => url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1'
#    91.0 (21-Aug-2013 to 04-Apr-2014)
#    => url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0'
#    90.9 (3-Jan-2011 to 20-Aug-2013)
#    => url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9'
#    90.8 (7-May-2009 to 2-Jan-2011)
#    => url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8'
#    90.6 (18-Sep-2008 to 6-May-2009)
#    => url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6'
# Legacy Experiments
#    74.2 (Jun-2007 to Oct-2008)
#    90.3 (Apr-2007 to Sep-2008)
#    90.2 (Jan-2007 to Apr-2007)
#    60.5 (Nov-2003 to Dec-2006)


class HYCOM(object):

    def __init__(self, dates):

        try:
            from pydap.client import open_url
        except ImportError:
            print('... you are in bad luck !\n'
                  '  ~>  you need to have the pydap library for '
                  'python 3 installed,\n'
                  '  ~>  along with its dependencies')
            sys.exit(0)
        # ~~~~ Initialisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.moddates = [datetime(*dates[0]), datetime(*dates[1])]
        self.slf2d = None
        self.slf3d = None
        self.hycomdata = None
        self.hycomilon = None
        self.hycomilat = None
        self.zplan = None

        # ~~~~ Time records ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Extract HYCOM time records\n')
        hycomurls = [
            'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2',
            'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1',
            'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0',
            'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9',
            'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8',
            'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6'
                    ]

        self.experiments = []
        for hycomurl in hycomurls:
            success = False
            attempt = 0
            while not success:
                try:
                    success = True
                    hycomdata = open_url(hycomurl)
                    nit = hycomdata['Date'].shape[0]
                    # /!\ the following print statement requires __future__
                    print('        x '+str(nit)+' records from '+hycomurl,
                          end='')
                    its = []
                    ats = []
                    z = zip(hycomdata['Date'][0:nit], range(nit))
                except Exception:
                    if attempt == 4:
                        raise TelemacException("Could not access hycom data\
                                                (Maybe proxy issue ?)")
                    success = False
                    attempt += 1
                    print(' ... re-attempting {}'.format(attempt))

            for hycomdate, itime in z:
                date = datetime(int(str(hycomdate)[0:4]),
                                int(str(hycomdate)[4:6]),
                                int(str(hycomdate)[6:8]))
                # /!\ the following print statement requires __future__
                if itime == 0:
                    print(' from: '+str(date), end='')
                if itime == nit-1:
                    print(' to: '+str(date))
                if self.moddates[0] <= date and date <= self.moddates[1]:
                    its.append(itime)
                    ats.append(date)
            if its != []:
                self.experiments.append((hycomdata, nit, its, ats, hycomurl))

        print('\n')

    def get_header_hycom(self, bounds):

        # ~~> inheritence
        self.slf3d = Selafin('')     # slf3d
        self.slf2d = Selafin('')     # slf2d surface

        print('     +> Set Selafin Variables')
        self.slf3d.title = ''
        self.slf3d.nbv1 = 6
        self.slf3d.nvar = 6
        self.slf3d.varindex = range(self.slf3d.nvar)
        self.slf3d.varnames = ['ELEVATION Z     ',
                               'SALINITY        ', 'TEMPERATURE     ',
                               'VELOCITY U      ', 'VELOCITY V      ',
                               'VELOCITY W      ']
        self.slf3d.varunits = ['M               ',
                               'G/L             ', 'DEGREES         ',
                               'M/S             ', 'M/S             ',
                               'M/S             ']
        self.slf2d.title = self.slf3d.title
        self.slf2d.nbv1 = self.slf3d.nbv1 + 1
        self.slf2d.nvar = self.slf3d.nvar + 1
        self.slf2d.varindex = range(self.slf2d.nvar)
        self.slf2d.varnames = self.slf3d.varnames[0:-1]
        self.slf2d.varnames.append('EMP             ')
        self.slf2d.varnames.append('QTOT            ')
        self.slf2d.varunits = self.slf3d.varunits[0:-1]
        self.slf2d.varunits.append('???             ')
        self.slf2d.varunits.append('???             ')
        # ~~> server access,
        #     get the grid and header from the latest experiment
        self.hycomdata = self.experiments[0][0]

        # ~~~~ Grid coordinates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        success = False
        while not success:
            try:
                success = True
                # ~~> the whole of the 2D grid sizes
                print('     +> Extract HYCOM sizes')
                nx1d = self.hycomdata['X'].shape[0]
                ny1d = self.hycomdata['Y'].shape[0]
                print('     +> Extract HYCOM mesh')
                lonx1d = self.hycomdata['Longitude']['Longitude']\
                    .data[0, 0:nx1d].ravel() % 360
                lat_y1d = self.hycomdata['Latitude']['Latitude']\
                    .data[0:ny1d, 0].ravel()
            except Exception:
                success = False
                print(' ... re-attempting ')
        # ~~> lat,lon correction
        for i in range(nx1d):
            if lonx1d[i] > 180:
                lonx1d[i] = lonx1d[i] - 360.0
        for i in range(2172, ny1d):
            lat_y1d[i] = 47.0 + (i-2172)/18.0
        # ~~> subset for the Selafin
        print('     +> Set Selafin mesh')
        self.hycomilon = \
            np.where((lonx1d >= bounds[0][1])*(lonx1d <= bounds[1][1]))[0]
        self.hycomilat = \
            np.where((lat_y1d >= bounds[0][0])*(lat_y1d <= bounds[1][0]))[0]
        x = lonx1d[self.hycomilon]
        y = lat_y1d[self.hycomilat]
        nx1d = len(x)
        ny1d = len(y)

        # ~~~~ MESH sizes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~> 3D
        success = False
        while not success:
            try:
                success = True
                print('     +> Set Selafin sizes')
                self.slf3d.nplan = self.hycomdata['Depth'].shape[0]
                # I do not know any other way
                self.zplan = self.hycomdata['Depth'][0:self.slf3d.nplan][::-1]
            except Exception:
                success = False
                print(' ... re-attempting ')
        self.slf3d.ndp2 = 3
        self.slf3d.ndp3 = 6
        self.slf3d.npoin2 = nx1d * ny1d
        self.slf3d.npoin3 = self.slf3d.npoin2 * self.slf3d.nplan
        self.slf3d.nelem2 = 2*(nx1d-1)*(ny1d-1)
        self.slf3d.nelem3 = self.slf3d.nelem2 * (self.slf3d.nplan-1)
        self.slf3d.iparam = [0, 0, 0, 0, 0, 0, self.slf3d.nplan, 0, 0, 0]
        # ~~> 2D
        self.slf2d.nplan = 1
        self.slf2d.ndp2 = self.slf3d.ndp2
        self.slf2d.ndp3 = self.slf2d.ndp2
        self.slf2d.npoin2 = self.slf3d.npoin2
        self.slf2d.npoin3 = self.slf2d.npoin2
        self.slf2d.nelem2 = self.slf3d.nelem2
        self.slf2d.nelem3 = self.slf2d.nelem2
        self.slf2d.iparam = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

        print('     +> Set Selafin mesh')
        self.slf3d.meshx = np.tile(x, ny1d).reshape(ny1d, nx1d).itime.ravel()
        self.slf3d.meshy = np.tile(y, nx1d)
        self.slf2d.meshx = self.slf3d.meshx[0:self.slf2d.npoin2]
        self.slf2d.meshy = self.slf3d.meshy[0:self.slf2d.npoin2]

        # ~~~~ Connectivity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Set Selafin ikle')
        ielem = 0
        pbar = ProgressBar(maxval=self.slf3d.nelem3).start()
        self.slf3d.ikle3 = np.zeros((self.slf3d.nelem3, self.slf3d.ndp3),
                                    dtype=np.int)
        for k in range(1, self.slf3d.nplan):
            for i in range(1, nx1d):
                for j in range(1, ny1d):
                    ipoin = (i-1)*ny1d + j - 1 + (k-1)*self.slf3d.npoin2
                    # ~~> first prism
                    self.slf3d.ikle3[ielem][0] = ipoin
                    self.slf3d.ikle3[ielem][1] = ipoin + ny1d
                    self.slf3d.ikle3[ielem][2] = ipoin + 1
                    self.slf3d.ikle3[ielem][3] = ipoin + self.slf3d.npoin2
                    self.slf3d.ikle3[ielem][4] = ipoin + ny1d + \
                        self.slf3d.npoin2
                    self.slf3d.ikle3[ielem][5] = ipoin + 1 + self.slf3d.npoin2
                    ielem = ielem + 1
                    pbar.update(ielem)
                    # ~~> second prism
                    self.slf3d.ikle3[ielem][0] = ipoin + ny1d
                    self.slf3d.ikle3[ielem][1] = ipoin + ny1d + 1
                    self.slf3d.ikle3[ielem][2] = ipoin + 1
                    self.slf3d.ikle3[ielem][3] = ipoin + ny1d + \
                        self.slf3d.npoin2
                    self.slf3d.ikle3[ielem][4] = ipoin + ny1d + 1 + \
                        self.slf3d.npoin2
                    self.slf3d.ikle3[ielem][5] = ipoin + 1 + self.slf3d.npoin2
                    ielem = ielem + 1
                    pbar.update(ielem)
        pbar.finish()
        self.slf2d.ikle3 = np.compress(np.repeat([True, False],
                                                 self.slf2d.ndp2),
                                       self.slf3d.ikle3[0:self.slf3d.nelem2],
                                       axis=1)

        # ~~~~ Boundaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Set Selafin ipobO')
        pbar = ProgressBar(maxval=nx1d+ny1d).start()
        self.slf3d.ipob3 = np.zeros(self.slf3d.npoin3, dtype=np.int)
        # ~~> along the x-axis (lon)
        for i in range(nx1d):
            for k in range(1, self.slf3d.nplan+1):
                ipoin = i*ny1d + (k-1)*(2*nx1d+2*ny1d-4)
                self.slf3d.ipob3[ipoin] = i + 1 + (k-1)*(2*nx1d+2*ny1d-4)
                ipoin = i*ny1d - 1 + (k-1)*(2*nx1d+2*ny1d-4)
                self.slf3d.ipob3[ipoin] = 2*nx1d+(ny1d-2) - i + \
                    (k-1)*(2*nx1d+2*ny1d-4)
            pbar.update(i)
        # ~~> along the y-axis (alt)
        for i in range(1, ny1d):
            for k in range(1, self.slf3d.nplan+1):
                ipoin = i + (k-1)*(2*nx1d+2*ny1d-4)
                self.slf3d.ipob3[ipoin] = 2*nx1d + 2*(ny1d-2) - i + 1 + \
                    (k-1)*(2*nx1d+2*ny1d-4)
                ipoin = ny1d*(nx1d-1) + i + (k-1)*(2*nx1d+2*ny1d-4)
                self.slf3d.ipob3[ipoin] = nx1d + i + (k-1)*(2*nx1d+2*ny1d-4)
            pbar.update(i+nx1d)
        pbar.finish()
        self.slf2d.ipob3 = self.slf3d.ipob3[0:self.slf3d.npoin2]

    def put_content(self, root_name, only_2d):

        nbar = 0
        for exp in self.experiments:
            nbar += len(exp[2])
        ilat = [self.hycomilat[0], self.hycomilat[-1]+1]
        ilon = [self.hycomilon[0], self.hycomilon[-1]+1]

        # ~~~~ Time records ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Extract HYCOM time records')
        self.slf3d.tags = {'times': []}
        # ~~~~ Start Date and Time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.slf3d.tags['times'] = 86400.0 * np.arange(nbar)
        self.slf2d.tags = {'times': self.slf3d.tags['times']}
        self.slf3d.datetime = self.experiments[-1][3][0].timetuple()[0:6]
        self.slf2d.datetime = self.slf3d.datetime
        self.slf3d.iparam[9] = 1
        self.slf2d.iparam[9] = 1

        print('     +> Write Selafin headers')
        if not only_2d:
            self.slf3d.fole = {}
            self.slf3d.fole.update({'hook': open('t3d_'+root_name, 'wb')})
            self.slf3d.fole.update({'name': 't3d_'+root_name})
            self.slf3d.fole.update({'endian': ">"})     # big endian
            self.slf3d.fole.update({'float': ('f', 4)})  # single precision
            self.slf3d.append_header_slf()
        self.slf2d.fole = {}
        self.slf2d.fole.update({'hook': open('t2d_'+root_name, 'wb')})
        self.slf2d.fole.update({'name': 't2d_'+root_name})
        self.slf2d.fole.update({'endian': ">"})     # big endian
        self.slf2d.fole.update({'float': ('f', 4)})  # single precision
        self.slf2d.append_header_slf()
        # ~~~~ Time loop(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        var3d = np.zeros(self.slf3d.npoin3, dtype=np.float)
        var2d = np.zeros(self.slf2d.npoin3, dtype=np.float)
        print('     +> Write Selafin cores')
        ibar = 0
        pbar = ProgressBar(maxval=10*nbar).start()
        for exp in self.experiments[::-1]:
            hycomdata = exp[0]
            i_1 = min(exp[2])
            i_2 = max(exp[2])+1
            for itime in range(i_1, i_2):
                # ~~> time stamp
                pbar.write('        x '+str(exp[3][itime-i_1]), 10*ibar+0)
                pbar.update(10*ibar+0)
                if not only_2d:
                    self.slf3d.append_core_time_slf(ibar)
                self.slf2d.append_core_time_slf(ibar)

                # ~~> HYCOM variable extraction
                #     ( 1L:times, 33L:layers, yyL:ny1d, xxL:nx1d )

                # ~~> ELEVATION
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - ssh', 10*ibar+1)
                        data = hycomdata['ssh']['ssh']
                        v2d = np.swapaxes(data.data[itime, ilat[0]:ilat[1],
                                                    ilon[0]:ilon[1]][0],
                                          0, 1).ravel()
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                    I failed ...',
                                   10*ibar)
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    var3d = - np.tile(self.zplan, self.slf3d.npoin2)\
                                .reshape(self.slf3d.npoin2, self.slf3d.nplan)\
                                .itime.ravel()
                    var3d[self.slf3d.npoin3-self.slf3d.npoin2:] = var2d
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.update(10*ibar+1)

                # ~~> SALINITY
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - surface_salinity_trend',
                                   10*ibar+2)
                        data = hycomdata['surface_salinity_trend'][
                                         'surface_salinity_trend']
                        v2d = np.swapaxes(
                                  data.data[itime, ilat[0]:ilat[1],
                                            ilon[0]:ilon[1]][0],
                                  0, 1).ravel()
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                            I failed ...',
                                   10*ibar)
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                pbar.update(10*ibar+2)
                if not only_2d:
                    success = False
                    while not success:
                        try:
                            success = True
                            pbar.write('             - salinity', 10*ibar+3)
                            data = hycomdata['salinity']['salinity']
                            var = np.swapaxes(
                                    data.data[itime, 0:self.slf3d.nplan,
                                              ilat[0]:ilat[1],
                                              ilon[0]:ilon[1]][0],
                                    1, 2)
                        except Exception:
                            success = False
                            pbar.write(
                                  '    ... re-attempting because I failed ...',
                                  10*ibar)
                    v3d = var[::-1].ravel()
                    var3d = np.where(v3d < 10000, v3d, 0.0)
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.update(10*ibar+3)

                # ~~> TEMPERATURE
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - surface_temperature_trend',
                                   10*ibar+4)
                        data = hycomdata['surface_temperature_trend'][
                                         'surface_temperature_trend']
                        v2d = np.swapaxes(
                                  data.data[itime, ilat[0]:ilat[1],
                                            ilon[0]:ilon[1]][0],
                                  0, 1).ravel()
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                            I failed ...',
                                   10*ibar)
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                pbar.update(10*ibar+4)
                if not only_2d:
                    success = False
                    while not success:
                        try:
                            success = True
                            pbar.write('             - temperature', 10*ibar+5)
                            data = hycomdata['temperature']['temperature']
                            var = np.swapaxes(
                                      data.data[itime, 0:self.slf3d.nplan,
                                                ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0],
                                      1, 2)
                        except Exception:
                            success = False
                            pbar.write(
                                 '    ... re-attempting because I failed ...',
                                 10*ibar)
                    v3d = var[::-1].ravel()
                    var3d = np.where(v3d < 10000, v3d, 0.0)
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.update(10*ibar+5)

                # ~~> VELOCITY U
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - u-velocity', 10*ibar+6)
                        data = hycomdata['u']['u']
                        if only_2d:
                            var = np.swapaxes(data.data[itime, 0:1,
                                                        ilat[0]:ilat[1],
                                                        ilon[0]:ilon[1]][0],
                                              1, 2)
                        else:
                            var = np.swapaxes(data.data[itime,
                                                        0:self.slf3d.nplan,
                                                        ilat[0]:ilat[1],
                                                        ilon[0]:ilon[1]][0],
                                              1, 2)
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                            I failed ...',
                                   10*ibar)
                v2d = var[0].ravel()
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    v3d = var[::-1].ravel()
                    var3d = np.where(v3d < 10000, v3d, 0.0)
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.update(10*ibar+6)

                # ~~> VELOCITY V
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - v-velocity', 10*ibar+7)
                        data = hycomdata['v']['v']
                        if only_2d:
                            var = np.swapaxes(data.data[itime, 0:1,
                                                        ilat[0]:ilat[1],
                                                        ilon[0]:ilon[1]][0],
                                              1, 2)
                        else:
                            var = np.swapaxes(data.data[itime,
                                                        0:self.slf3d.nplan,
                                                        ilat[0]:ilat[1],
                                                        ilon[0]:ilon[1]][0],
                                              1, 2)
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                            I failed ...',
                                   10*ibar)
                v2d = var[0].ravel()
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    v3d = var[::-1].ravel()
                    var3d = np.where(v3d < 10000, v3d, 0.0)
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.update(10*ibar+7)

                # ~~> VELOCITY W
                if not only_2d:
                    var3d = 0. * var3d
                    self.slf3d.append_core_vars_slf([var3d])

                # ~~> EMP ???
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - emp', 10*ibar+8)
                        data = hycomdata['emp']['emp']
                        v2d = np.swapaxes(data.data[itime, ilat[0]:ilat[1],
                                                    ilon[0]:ilon[1]][0],
                                          0, 1).ravel()
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                            I failed ...',
                                   10*ibar)
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                pbar.update(10*ibar+8)

                # ~~> TEMPERATURE
                success = False
                while not success:
                    try:
                        success = True
                        pbar.write('             - qtot',
                                   10*ibar+9)
                        data = hycomdata['qtot']['qtot']
                        v2d = np.swapaxes(data.data[itime, ilat[0]:ilat[1],
                                                    ilon[0]:ilon[1]][0],
                                          0, 1).ravel()
                    except Exception:
                        success = False
                        pbar.write('    ... re-attempting because\
                                            I failed ...',
                                   10*ibar)
                var2d = np.where(v2d < 10000, v2d, 0.0)
                self.slf2d.append_core_vars_slf([var2d])
                pbar.update(10*ibar+9)

                ibar += 1

        pbar.finish()
        if not only_2d:
            self.slf3d.fole['hook'].close()
        self.slf2d.fole['hook'].close()

    def __del__(self):
        pass


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
def hycom2srf_parser(subparser):
    """
    Adding options for hycom2srf

    @param subparser (ArgumentParser) The argument parser to update

    @return (ArgumentParser) the updated argument parser
    """
    parser = subparser.add_parser('hycom2srf',
                                  help='Extracting HYCOM data into\
                                        a serafin-format file')
    parser.add_argument(
        "-r", "--root", dest="root_name", default='hycom.slf', required=True,
        help="root name used for the output")
    parser.add_argument(
        "-f", "--from", dest="tfrom", default=None, required=True,
        help="specify the first date included (1972-13-07)")
    parser.add_argument(
        "-s", "--stop", dest="tstop", default=None, required=True,
        help="specify the last date included (1980-12-31)")
    parser.add_argument(
        "--bl", dest="blcorner", default=None, required=True,
        help="specify the bottom left corner (25,-117)")
    parser.add_argument(
        "--tr", dest="trcorner", default=None, required=True,
        help="specify the top right corner (27,-110)")
    parser.add_argument(
        "--2d", action="store_true", dest="t2d", default=False,
        help="if there, produces on the 2D file")

    return subparser


def hycom2srf(tfrom, tstop, blcorner, trcorner, root_name, only_2d):
    """
    Convertion form HYCOM data to Serafin format
    """
    # Arbitrary 6-day period
    period = [[], []]
    if tfrom is not None:
        for i in tfrom.split('-'):
            period[0].append(int(i))
    else:
        raise TelemacException(
                '... could not find your from date. '
                'Please use --from option '
                '(- delimited, no spaces).\n\n')
    if tstop is not None:
        for i in tstop.split('-'):
            period[1].append(int(i))
    else:
        raise TelemacException(
                '... could not find your stop date. '
                'Please use --stop option '
                '(- delimited, no spaces).\n\n')

    # arbitrary box (small pieve of the atlantic side of Mexico)
    modelbox = [[], []]
    if blcorner is not None:
        for i in blcorner.split(','):
            modelbox[0].append(int(i))
    else:
        raise TelemacException(
             '... could not find your bounding box bottom left corner. '
             'Please use --bl option (, delimited, no spaces).\n\n')
    if trcorner is not None:
        for i in trcorner.split(','):
            modelbox[1].append(int(i))
    else:
        raise TelemacException(
            '... could not find your bounding box top right corner. '
            'Please use --tr option (, delimited, no spaces).\n\n')

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Convert to Selafin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('\n\n'+'~'*72+'\n')
    print('\nProcessing header (mesh, connectivity, etc.)\n')
    hy2slf = HYCOM(period)
    hy2slf.get_header_hycom(modelbox)

    print('\n\n'+'~'*72+'\n')
    print('\nProcessing core variables (time record, variables, etc.)\n')
    tic = time.time()
    print('\nExtraction start time:   ' +
          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(tic)))
    hy2slf.put_content(root_name, only_2d)
    toc = time.time()
    print('\nExtraction end time:     ' +
          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(toc)))
    print('___________\nDuration:     '+str(int(toc-tic))+' seconds\n')


__author__ = "Juliette C.E. Parisi"
__date__ = "$13-Nov-2013 08:51:29$"
