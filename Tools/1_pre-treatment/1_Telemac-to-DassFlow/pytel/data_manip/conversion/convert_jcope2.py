#!/usr/bin/env python
r"""@author Sebastien E. Bourban and Noemie Durand

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
    Tools for handling conversions to-from JCOPE2 server files

    @details
    Contains server read functions to convert to Selafin file

    @history 30/04/2013 -- Sebastien E. Bourban
    Frist draft, also the final version.

    @history 11/09/2013 -- Noemie Durand
    Modification of the 3D output file, so the elevation of the top
      layer is the same as the elevation in the 2D file.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import time
from datetime import datetime, timedelta
import numpy as np
# ~~> dependencies towards other pytel scripts
# ~~> dependencies towards other modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
# ~~> dependencies towards other modules


# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#
# Current Experiments
#    => url = 'http://apdrc.soest.hawaii.edu/dods/public_data/FRA-JCOPE2'
# Other Experiments
#    => url = 'http://apdrc.soest.hawaii.edu/dapper/public_data/
# topography/io_bathy/modified_etopo2.nc'

class Jcope2(object):

    def __init__(self, dates):

        try:
            from pydap.client import open_url
        except Exception as excpt:
            raise TelemacException('... you are in bad luck !\n'
                                   '  ~>  you need the pydap library\
                                    unzipped locally\n'
                                   '{}'.format(excpt))
        # ~~~~ Initialisation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.moddates = [datetime(*dates[0]), datetime(*dates[1])]
        jcope2vars = ['el', 'itime', 's', 'u', 'v']
        # /!\ unknown convertion of time records into dates
        jcope2date = [1993, 1, 1]
        jcope2root = \
            'http://apdrc.soest.hawaii.edu/dods/public_data/FRA-JCOPE2'
        self.slf2d = None
        self.slf3d = None
        self.jcope2ilon = None
        self.jcope2ilat = None
        self.zplan = None
        self.mask2 = None
        self.mask3 = None

        # ~~~~ Time records ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Extract JCOPE2 time records\n')
        self.experiments = []
        experiment = {}  # /!\ only one time period covered at this stage
        for jvar in jcope2vars:
            jcope2url = jcope2root+'/'+jvar
            jcope2data = open_url(jcope2url)
            nit = jcope2data['time'].shape[0]
            # /!\ the following print statement requires __future__
            print('        x '+str(nit)+' records from '+jcope2url, end='')
            its = []
            ats = []
            for itime in range(nit):
                date = datetime(jcope2date[0], jcope2date[1], jcope2date[2]) +\
                     timedelta(itime)
                if itime == 0:
                    print(' from: ' + str(date), end='')
                if itime == nit-1:
                    print(' to: ' + str(date))
                if self.moddates[0] <= date and date <= self.moddates[1]:
                    its.append(itime)
                    ats.append(date)
            if its != []:
                for ivar in list(jcope2data.keys()):
                    if ivar not in ['time', 'lev', 'lat', 'lon']:
                        experiment.update({ivar: jcope2data[ivar]})
            else:
                raise TelemacException(
                    '... I could not find the time to do your work\n'
                    '  ~>  you may need to select a different time period\n')
        self.experiments.append((experiment, nit, its, ats))
        print('\n')

    def get_header_jcope2(self, bounds):

        # ~~> inheritence
        self.slf3d = Selafin('')     # slf3d
        self.slf2d = Selafin('')     # slf2d surface

        print('     +> Set Selafin Variables')
        self.slf3d.title = ''
        self.slf3d.nbv1 = 6
        self.slf3d.nvar = 6
        self.slf3d.varindex = list(range(self.slf3d.nvar))
        self.slf3d.varnames = ['ELEVATION Z     ',
                               'SALINITY        ', 'TEMPERATURE     ',
                               'VELOCITY U      ', 'VELOCITY V      ',
                               'VELOCITY W      ']
        self.slf3d.varunits = ['M               ',
                               '                ', '                ',
                               'M/S             ', 'M/S             ',
                               'M/S             ']
        self.slf2d.title = self.slf3d.title
        self.slf2d.nbv1 = self.slf3d.nbv1 - 1
        self.slf2d.nvar = self.slf2d.nbv1
        self.slf2d.varindex = list(range(self.slf2d.nvar))
        self.slf2d.varnames = self.slf3d.varnames[0:-1]
        self.slf2d.varunits = self.slf3d.varunits[0:-1]

        # ~~~~ Grid coordinates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # ~~> the whole of the 2D grid sizes
        print('     +> Extract JCOPE2 sizes')
        # /!\ 'itime' gives me access to nplan in 3D
        jcope2data = self.experiments[0][0]['temp']
        nx1d = jcope2data['lon'].shape[0]
        ny1d = jcope2data['lat'].shape[0]
        print('     +> Extract JCOPE2 mesh')
        lon_x1d = jcope2data['lon'].data[0:nx1d]
        lat_y1d = jcope2data['lat'].data[0:ny1d]
        # ~~> no correction for lat,lon
        # ~~> subset for the Selafin
        print('     +> Set Selafin mesh')
        self.jcope2ilon = \
            np.where((lon_x1d >= bounds[0][1]) * (lon_x1d <= bounds[1][1]))[0]
        self.jcope2ilat = \
            np.where((lat_y1d >= bounds[0][0]) * (lat_y1d <= bounds[1][0]))[0]
        x = lon_x1d[self.jcope2ilon]
        y = lat_y1d[self.jcope2ilat]
        nx1d = len(x)
        ny1d = len(y)

        # ~~~~ MESH sizes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Set Selafin sizes')
        # ~~> 3D
        self.slf3d.nplan = jcope2data['lev'].shape[0]
        # I do not know any other way
        self.zplan = jcope2data['lev'][0:self.slf3d.nplan][::-1]
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

        # ~~~~ Connectivity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Set the default Selafin ikle 3D')
        ielem = 0
        pbar = ProgressBar(maxval=self.slf3d.nelem3).start()
        self.slf3d.ikle3 = np.zeros((self.slf3d.nelem3, self.slf3d.ndp3),
                                    dtype=np.int)
        npoin2 = self.slf3d.npoin2
        for k in range(1, self.slf3d.nplan):
            for i in range(1, nx1d):
                for j in range(1, ny1d):
                    ipoin = (i-1)*ny1d + j - 1 + (k-1)*npoin2
                    # ~~> first prism
                    self.slf3d.ikle3[ielem][0] = ipoin
                    self.slf3d.ikle3[ielem][1] = ipoin + ny1d
                    self.slf3d.ikle3[ielem][2] = ipoin + 1
                    self.slf3d.ikle3[ielem][3] = ipoin + npoin2
                    self.slf3d.ikle3[ielem][4] = ipoin + ny1d + npoin2
                    self.slf3d.ikle3[ielem][5] = ipoin + 1 + npoin2
                    ielem = ielem + 1
                    pbar.update(ielem)
                    # ~~> second prism
                    self.slf3d.ikle3[ielem][0] = ipoin + ny1d
                    self.slf3d.ikle3[ielem][1] = ipoin + ny1d + 1
                    self.slf3d.ikle3[ielem][2] = ipoin + 1
                    self.slf3d.ikle3[ielem][3] = ipoin + ny1d + npoin2
                    self.slf3d.ikle3[ielem][4] = ipoin + ny1d + 1 + npoin2
                    self.slf3d.ikle3[ielem][5] = ipoin + 1 + npoin2
                    ielem = ielem + 1
                    pbar.update(ielem)
        pbar.finish()
        self.slf2d.ikle3 = np.compress([True, True, True, False, False, False],
                                       self.slf3d.ikle3[0:self.slf3d.nelem2],
                                       axis=1)

        # ~~~~ Boundaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Set Selafin ipobO')
        pbar = ProgressBar(maxval=nx1d+ny1d).start()
        ipob2 = np.zeros(self.slf3d.npoin2, dtype=np.int)
        # ~~> along the x-axis (lon)
        for i in range(nx1d):
            ipoin = i*ny1d
            ipob2[ipoin] = i + 1
            ipoin = i*ny1d - 1
            ipob2[ipoin] = 2*nx1d+(ny1d-2) - i
            pbar.update(i)
        # ~~> along the y-axis (alt)
        for i in range(1, ny1d):
            ipoin = i
            ipob2[ipoin] = 2*nx1d + 2*(ny1d-2) - i + 1
            ipoin = ny1d*(nx1d-1) + i
            ipob2[ipoin] = nx1d + i
            pbar.update(i+nx1d)
        pbar.finish()

        # ~~~~ Connectivity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # /!\ 'el' gives me access to the real mesh removing
        #       elements with -99 values
        print('     +> Mask the non-values from the Selafin ikle')
        jcope2data = self.experiments[0][0]['el']
        var = np.swapaxes(jcope2data['el']
                          .data[0, 0, self.jcope2ilat[0]:self.jcope2ilat[-1]+1,
                                self.jcope2ilon[0]:self.jcope2ilon[-1]+1][0],
                          1, 2).ravel()
        # ~> the elements you wish to keep
        array = np.compress(var > -99, np.arange(len(var)))
        array_1d = np.in1d(self.slf2d.ikle3, array)
        array_sum = np.sum(array_1d.reshape(self.slf2d.nelem2,
                                            self.slf2d.ndp2),
                           axis=1)
        mask2 = self.slf2d.ikle3[np.where(array_sum == 3)]

        self.slf2d.nelem2 = len(mask2)
        self.slf2d.nelem3 = self.slf2d.nelem2
        self.slf3d.nelem2 = self.slf2d.nelem2
        self.slf3d.nelem3 = self.slf3d.nelem2 * (self.slf3d.nplan-1)

        # ~~> re-numbering ikle2 as a local connectivity matrix
        knolg, _ = np.unique(np.ravel(mask2), return_index=True)
        knogl = dict(list(zip(knolg, list(range(len(knolg))))))
        self.mask2 = np.in1d(np.arange(len(var)), knolg)
        self.mask3 = np.tile(self.mask2, self.slf3d.nplan)
        self.slf2d.ikle2 = - np.ones_like(mask2, dtype=np.int)
        for k in range(len(mask2)):
            self.slf2d.ikle2[k] = [knogl[mask2[k][0]],
                                   knogl[mask2[k][1]],
                                   knogl[mask2[k][2]]]

        self.slf3d.npoin2 = len(knolg)
        self.slf3d.npoin3 = self.slf3d.npoin2 * self.slf3d.nplan
        self.slf2d.npoin2 = self.slf3d.npoin2
        self.slf2d.npoin3 = self.slf2d.npoin2

        # ~~> re-connecting the upper floors
        self.slf2d.ikle3 = self.slf2d.ikle2
        self.slf3d.ikle2 = self.slf2d.ikle2

        init = self.slf2d.npoin2*np.arange(self.slf3d.nplan-1)
        size = self.slf2d.nelem2*self.slf3d.ndp3
        self.slf3d.ikle3 = \
            np.repeat(init, size).reshape((self.slf2d.nelem2 *
                                          (self.slf3d.nplan-1),
                                           self.slf3d.ndp3)) + \
            np.tile(np.add(np.tile(self.slf2d.ikle2, 2),
                           np.repeat(self.slf2d.npoin2*np.arange(2),
                                     self.slf3d.ndp2)),
                          (self.slf3d.nplan-1, 1))

        # ~~~~ Boundaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.slf2d.ipob2 = ipob2[self.mask2]
        self.slf2d.ipob3 = self.slf2d.ipob2
        self.slf3d.ipob2 = self.slf2d.ipob2
        self.slf3d.ipob3 = np.ravel(np.add(
                  np.repeat(self.slf2d.ipob2, self.slf3d.nplan)
                    .reshape((self.slf2d.npoin2, self.slf3d.nplan)),
                  self.slf2d.npoin2*np.arange(self.slf3d.nplan)).T)

        # ~~~~ Mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Set Selafin mesh')
        self.slf3d.meshx = np.tile(x, ny1d).reshape(ny1d, nx1d)\
            .T.ravel()[self.mask2] + 0.042
        self.slf3d.meshy = np.tile(y, nx1d)[self.mask2] + 0.042
        self.slf2d.meshx = self.slf3d.meshx
        self.slf2d.meshy = self.slf3d.meshy

    def put_content(self, root_name, only_2d):

        nbar = 0
        for exp in self.experiments:
            nbar += len(exp[2])
        ilat = [self.jcope2ilat[0], self.jcope2ilat[-1]+1]
        ilon = [self.jcope2ilon[0], self.jcope2ilon[-1]+1]

        # ~~~~ Time records ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('     +> Extract JCOPE2 time records')
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
        print('     +> Write Selafin cores')
        ibar = 0
        pbar = ProgressBar(maxval=6*nbar).start()
        for exp in self.experiments:
            jcope2data = exp[0]
            i_1 = min(exp[2])
            i_2 = max(exp[2])+1
            for itime in range(i_1, i_2):
                # ~~> time stamp
                pbar.write('        x '+str(exp[3][itime-i_1]), 6*ibar+0)
                pbar.update(6*ibar+0)
                if not only_2d:
                    self.slf3d.append_core_time_slf(ibar)
                self.slf2d.append_core_time_slf(ibar)

                # ~~> ELEVATION
                data = jcope2data['el']['el']
                var2d = np.swapaxes(data.data[itime, 0, ilat[0]:ilat[1],
                                              ilon[0]:ilon[1]][0],
                                    1, 2).ravel()[self.mask2]
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    var3d = - np.tile(self.zplan, self.slf3d.npoin2)\
                                .reshape(self.slf3d.npoin2, self.slf3d.nplan)\
                                .T.ravel()
                    var3d[self.slf3d.npoin3-self.slf3d.npoin2:] = var2d
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.write('             - elevation', 6*ibar+1)
                pbar.update(6*ibar+1)

                # ~~> SALINITY
                data = jcope2data['salt']['salt']
                if only_2d:
                    var = np.swapaxes(data.data[itime, 0:1, ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                else:
                    var = np.swapaxes(data.data[itime, 0:self.slf3d.nplan,
                                                ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0],
                                      1, 2)
                var2d = var[0].ravel()[self.mask2]
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    var3d = var[::-1].ravel()[self.mask3]
                    for ipoin in range(self.slf3d.npoin2):
                        for iplan in range(self.slf3d.nplan-1, 0, -1):
                            if var3d[ipoin+(iplan-1)*self.slf3d.npoin2] \
                                    < -99.0:
                                var3d[ipoin+(iplan-1)*self.slf3d.npoin2] = \
                                          var3d[ipoin+iplan*self.slf3d.npoin2]
                    self.slf3d.append_core_vars_slf([var3d])
                    pbar.write('             - salinity', 6*ibar+2)
                pbar.update(6*ibar+2)

                # ~~> TEMPERATURE
                data = jcope2data['temp']['temp']
                if only_2d:
                    var = np.swapaxes(data.data[itime, 0:1, ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                else:
                    var = np.swapaxes(data.data[itime, 0:self.slf3d.nplan,
                                                ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                var2d = var[0].ravel()[self.mask2]
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    var3d = var[::-1].ravel()[self.mask3]
                    for ipoin in range(self.slf3d.npoin2):
                        for iplan in range(self.slf3d.nplan-1, 0, -1):
                            if var3d[ipoin+(iplan-1)*self.slf3d.npoin2] \
                                    < -99.0:
                                var3d[ipoin+(iplan-1)*self.slf3d.npoin2] = \
                                          var3d[ipoin+iplan*self.slf3d.npoin2]
                    self.slf3d.append_core_vars_slf([var3d])
                    pbar.write('             - temperature', 6*ibar+3)
                pbar.update(6*ibar+3)

                # ~~> VELOCITY U
                data = jcope2data['u']['u']
                if only_2d:
                    var = np.swapaxes(data.data[itime, 0:1, ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                else:
                    var = np.swapaxes(data.data[itime, 0:self.slf3d.nplan,
                                                ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                var2d = var[0].ravel()[self.mask2]
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    var3d = var[::-1].ravel()[self.mask3]
                    for ipoin in range(self.slf3d.npoin2):
                        for iplan in range(self.slf3d.nplan-1, 0, -1):
                            if var3d[ipoin+(iplan-1)*self.slf3d.npoin2] \
                                    < -99.0:
                                var3d[ipoin+(iplan-1)*self.slf3d.npoin2] = \
                                          var3d[ipoin+iplan*self.slf3d.npoin2]
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.write('             - u-velocity', 6*ibar+4)
                pbar.update(6*ibar+4)

                # ~~> VELOCITY V
                data = jcope2data['v']['v']
                if only_2d:
                    var = np.swapaxes(data.data[itime, 0:1, ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                else:
                    var = np.swapaxes(data.data[itime, 0:self.slf3d.nplan,
                                                ilat[0]:ilat[1],
                                                ilon[0]:ilon[1]][0], 1, 2)
                var2d = var[0].ravel()[self.mask2]
                self.slf2d.append_core_vars_slf([var2d])
                if not only_2d:
                    var3d = var[::-1].ravel()[self.mask3]
                    for ipoin in range(self.slf3d.npoin2):
                        for iplan in range(self.slf3d.nplan-1, 0, -1):
                            if var3d[ipoin+(iplan-1)*self.slf3d.npoin2] \
                                    < -99.0:
                                var3d[ipoin+(iplan-1)*self.slf3d.npoin2] = \
                                          var3d[ipoin+iplan*self.slf3d.npoin2]
                    self.slf3d.append_core_vars_slf([var3d])
                pbar.write('             - v-velocity', 6*ibar+5)
                pbar.update(6*ibar+5)

                # ~~> VELOCITY W
                if not only_2d:
                    var3d = 0. * var3d
                    self.slf3d.append_core_vars_slf([var3d])

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

__author__ = "Sebastien E. Bourban"
__date__ = "$13-July-2014 08:51:29$"


def main():
    """ Main function of convertJCOPE2 """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+72*'~'+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
Download JCOPE2 data into a Selafin file
        '''),
        usage=' (--help for help)\n---------\n       =>  '
              '%(prog)s [options]\n---------',
        epilog=('''\nexamples:\n---------
1: Extract about 80 days of 2D results only from January 2, 1993.
        => convertJCOPE2.py --from 1993-01-02 --stop 1993-03-25 --bl 34,140
         --tr 41,147 -r jcope2-80d.slf --2d
2: Extract 4 months from both 2D and 3D dataset
        => convertJCOPE2.py --from 1993-01-02 --stop 1993-05-02 --bl 34,140
         --tr 41,147 -r jcope2-4m1993.slf
        => convertJCOPE2.py --from 1993-01-02 --stop 1993-05-02 --bl 10,108
         --tr 62,180 -r visu-jcope2-4m1993.slf
---------'''))
    parser.add_argument(
        "-r", "--root", dest="root_name", default='jcope2.slf', required=True,
        help="root name used for the output")
    parser.add_argument(
        "-f", "--from", dest="tfrom", default=None, required=True,
        help="specify the first date included (1972-13-07)")
    parser.add_argument(
        "-s", "--stop", dest="tstop", default=None, required=True,
        help="specify the last date included (1980-12-31)")
    parser.add_argument(
        "--bl", dest="blcorner", default=None, required=True,
        help="specify the bottom left corner (30,130)")
    parser.add_argument(
        "--tr", dest="trcorner", default=None, required=True,
        help="specify the top right corner (64,170)")
    parser.add_argument(
        "--2d", action="store_true", dest="t2d", default=False,
        help="if there, produces on the 2D file")
    options = parser.parse_args()

    # Arbitrary 6-day period
    period = [[], []]
    if options.tfrom is not None:
        for i in options.tfrom.split('-'):
            period[0].append(int(i))
    else:
        raise TelemacException('... could not find your from date. '
                               'Please use --from option '
                               '(- delimited, no spaces).\n\n')
    if options.tstop is not None:
        for i in options.tstop.split('-'):
            period[1].append(int(i))
    else:
        raise TelemacException('... could not find your stop date. '
                               'Please use --stop option '
                               '(- delimited, no spaces).\n\n')

    # arbitrary box (small pieve of the atlantic side of Mexico)
    modelbox = [[], []]
    if options.blcorner is not None:
        for i in options.blcorner.split(','):
            modelbox[0].append(int(i))
    else:
        raise TelemacException('... could not find your bounding box'
                               ' bottom left corner. Please use --bl option'
                               ' (, delimited, no spaces).\n\n')
    if options.trcorner is not None:
        for i in options.trcorner.split(','):
            modelbox[1].append(int(i))
    else:
        raise TelemacException('... could not find your bounding box'
                               ' top right corner. Please use --tr option '
                               '(, delimited, no spaces).\n\n')
    # root_name
    root_name = options.root_name

    # 2D file only
    only_2d = options.t2d

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Convert to Selafin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('\n\n'+'~'*72+'\n')
    print('\nProcessing header (mesh, connectivity, etc.)\n')
    jc2slf = Jcope2(period)
    jc2slf.get_header_jcope2(modelbox)

    print('\n\n'+'~'*72+'\n')
    print('\nProcessing core variables (time record, variables, etc.)\n')
    tic = time.time()
    print('\nExtraction start time:   ' +
          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(tic)))
    jc2slf.put_content(root_name, only_2d)
    toc = time.time()
    print('\nExtraction end time:     ' +
          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(toc)))
    print('___________\nDuration:     '+str(int(toc-tic))+' seconds\n')

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
