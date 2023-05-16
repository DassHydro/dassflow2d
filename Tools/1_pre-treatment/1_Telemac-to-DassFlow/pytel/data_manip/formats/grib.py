"""
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import numpy as np
from datetime import datetime, timedelta
import math
# ~~> dependencies towards other modules
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
#
try:
    import pygrib
    IMPORT_GRIB = True
except ImportError:
    IMPORT_GRIB = False


class Grib(object):

    def __init__(self, dataset, request, stream):

        self.request = request
        # ~~> inheritence
        self.slf2d = Selafin('')
        self.slf2d.title = ''
        self.slf2d.fole = {}
        self.dataset = []
        self.variables = []
        self.byrowdown = False
        self.maskx = None
        self.masky = None
        self.nx1d = None
        self.ny1d = None
        self.nb_direct = None
        self.nb_freq = None
        self.freq = None
        self.dirc = None

        print('   +> identifying relevant files, by variables')
        # ~~> filter requested variables
        self.variables = []
        found_dataset = []
        ibar = 0
        pbar = ProgressBar(maxval=len(dataset)).start()
        for data in dataset:
            grbs = pygrib.open(data)
            for grb in grbs:
                if str(grb.indicatorOfParameter) in request['param']\
                                                        .split('/'):
                    if data not in found_dataset:
                        found_dataset.append(data)
                    if grb.indicatorOfParameter not in self.variables:
                        self.variables.append(grb.indicatorOfParameter)
                    else:
                        break
            grbs.close()
            ibar += 1
            pbar.update(ibar)
        pbar.finish()
        if self.variables == []:
            raise TelemacException('... could not find\
                the requested valiables.\n\n')

        print('   +> sorting out timeline')
        # ~~>  checking consistency of origin of date and time
        for data in found_dataset:
            grbs = pygrib.open(data)
            for grb in grbs:
                at0 = pygrib.julian_to_datetime(grb.julianDay)
                break
            break
        print('      - start date and time', at0)
        self.slf2d.datetime = [d for d in at0.timetuple()[0:6]]
        # ~~>  recording times from origin of date and time
        ats = []
        dts = []
        ibar = 0
        pbar = ProgressBar(maxval=len(found_dataset)).start()
        for data in found_dataset:
            grbs = pygrib.open(data)
            for grb in grbs:
                date = str(grb.validityDate)
                ats.append(datetime(int(date[:4]), int(date[4:6]),
                                    int(date[6:])) +
                           timedelta(seconds=int(grb.validityTime*36)))
                dts.append((ats[-1]-at0).total_seconds())
                break
            ibar += 1
            pbar.update(ibar)
        pbar.finish()
        print('      - finish date and time', ats[-1])
        # ~~>  checking if the list is sorted
        if not all(ats[i] < ats[i+1] for i in range(len(ats)-1)):
            raise TelemacException(
                    '... your dataset is not sorted. '
                    'Here is the time profile in seconds:\n'
                    '{}\n\n'.format(repr(dts)))
        # ~~> filter requested times
        udates = [datetime(*[int(a) for a in d.split('-')])
                  for d in request['date'].split('/to/')]
        self.slf2d.tags = {'times': []}
        udates[1] = udates[1]+timedelta(hours=24.)
        for i in range(len(ats)):
            if udates[0] <= ats[i] and ats[i] <= udates[1]:
                self.slf2d.tags['times'].append(dts[i])
                self.dataset.append(found_dataset[i])
        times = self.slf2d.tags['times']
        print('   +> actual timeline')
        print('      - start date and time  ', at0+timedelta(seconds=times[0]))
        print('      - finish date and time ',
              at0+timedelta(seconds=times[-1]))

        # ~> Other initialisations
        self.typ = stream

        # ~~> spatial sizes
        print('   +> checking out sizes')
        grbs = pygrib.open(self.dataset[0])
        for grb in grbs:
            self.missing_value = grb.missingValue
            self.scale_values_by = grb.scaleValuesBy
            self.offset = grb.offset
            break
        grbs.close()

    def open_grib(self, file_name):

        self.slf2d.fole.update({'hook': open(file_name, 'wb')})
        self.slf2d.fole.update({'name': file_name})
        self.slf2d.fole.update({'endian': ">"})     # big endian
        self.slf2d.fole.update({'float': ('f', 4)})  # single precision

    def close_grib(self):
        self.slf2d.fole['hook'].close()

    def set_geometry(self):

        # ~~> header
        self.byrowdown = False

        # ~~> 2D grid
        print('   +> set the mesh and connectivity')
        x_1, y_1, x_2, y_2 = self.request['area'].split('/')
        grbs = pygrib.open(self.dataset[0])
        for grb in grbs:
            y, x = grb.latlons()
            self.maskx = np.logical_and(float(x_1) <= x[0], x[0] <= float(x_2))
            l_x = x[0][self.maskx]
            if not np.any(self.maskx):
                self.maskx = np.logical_and(float(x_1) <= x[0]-360.,
                                            x[0]-360. <= float(x_2))
                l_x = x[0][self.maskx]-360.
            if not np.any(self.maskx):
                raise TelemacException(
                    '... your spatial range seems out of bound:\n       '
                    'you asked for [ {} - {}], while x is:\n       '
                    '{}\n\n'.format(x_1, x_2, repr(x)))
            self.nx1d = len(l_x)
            self.masky = np.logical_and(float(y_1) <= y.T[0],
                                        y.T[0] <= float(y_2))
            l_y = y.T[0][self.masky]
            if not np.any(self.masky):
                raise TelemacException(
                    '... your spatial range seems out of bound:\n       '
                    'you asked for [ {} - {}], while x is:\n       '
                    '{}\n\n'.format(y_1, y_2, repr(y)))
            self.ny1d = len(l_y)
            if self.byrowdown:
                self.slf2d.meshx = np.ravel(np.tile(l_x, self.ny1d)
                                            .reshape(self.ny1d, self.nx1d))
                self.slf2d.meshy = np.ravel(np.tile(l_y, self.nx1d)
                                            .reshape(self.nx1d, self.ny1d).T)
            else:
                self.slf2d.meshx = np.ravel(np.tile(l_x, self.ny1d)
                                            .reshape(self.ny1d, self.nx1d).T)
                self.slf2d.meshy = np.ravel(np.tile(l_y, self.nx1d)
                                            .reshape(self.nx1d, self.ny1d))
            break
        grbs.close()

        self.slf2d.nplan = 1
        self.slf2d.ndp2 = 3
        self.slf2d.ndp3 = self.slf2d.ndp2
        self.slf2d.npoin2 = self.nx1d * self.ny1d
        self.slf2d.npoin3 = self.slf2d.npoin2
        self.slf2d.nelem2 = 2*(self.nx1d-1)*(self.ny1d-1)
        self.slf2d.nelem3 = self.slf2d.nelem2

        # ~~> Connectivity - numbered by rows
        ielem = 0
        pbar = ProgressBar(maxval=self.slf2d.nelem3).start()
        self.slf2d.ikle3 = np.zeros((self.slf2d.nelem3, self.slf2d.ndp3),
                                    dtype=np.int)
        if self.byrowdown:
            for j in range(1, self.ny1d):
                for i in range(1, self.nx1d):
                    ipoin = (j-1)*self.nx1d + i - 1
                    # ~~> first triangle
                    self.slf2d.ikle3[ielem][0] = ipoin
                    self.slf2d.ikle3[ielem][1] = ipoin + self.nx1d
                    self.slf2d.ikle3[ielem][2] = ipoin + 1
                    ielem = ielem + 1
                    pbar.update(ielem)
                    # ~~> second triangle
                    self.slf2d.ikle3[ielem][0] = ipoin + self.nx1d
                    self.slf2d.ikle3[ielem][1] = ipoin + self.nx1d + 1
                    self.slf2d.ikle3[ielem][2] = ipoin + 1
                    ielem = ielem + 1
                    pbar.update(ielem)
        else:
            for j in range(1, self.ny1d):
                for i in range(1, self.nx1d):
                    ipoin = j - 1 + (i - 1)*self.ny1d
                    # ~~> first triangle
                    self.slf2d.ikle3[ielem][0] = ipoin
                    self.slf2d.ikle3[ielem][1] = ipoin + 1
                    self.slf2d.ikle3[ielem][2] = ipoin + self.ny1d
                    ielem = ielem + 1
                    pbar.update(ielem)
                    # ~~> second triangle
                    self.slf2d.ikle3[ielem][0] = ipoin + self.ny1d
                    self.slf2d.ikle3[ielem][1] = ipoin + 1
                    self.slf2d.ikle3[ielem][2] = ipoin + self.ny1d + 1
                    ielem = ielem + 1
                    pbar.update(ielem)
        pbar.finish()

        # ~~> Boundaries
        self.slf2d.ipob3 = np.zeros(self.slf2d.npoin3, dtype=np.int)

        if self.byrowdown:
            # ~~> around the box
            for i in range(self.nx1d):
                ipoin = i
                self.slf2d.ipob3[i] = ipoin
            for i in range(self.nx1d):
                ipoin = self.nx1d + self.ny1d - 2 + i
                self.slf2d.ipob3[self.nx1d*self.ny1d-i-1] = ipoin
            for j in range(1, self.ny1d-1):
                ipoin = j + self.nx1d - 1
                self.slf2d.ipob3[(j+1)*self.nx1d-1] = ipoin
            for j in range(1, self.ny1d-1):
                ipoin = self.ny1d + 2 * self.nx1d + j - 3
                self.slf2d.ipob3[self.nx1d*self.ny1d-j*self.nx1d-self.nx1d] = \
                    ipoin
        else:
            # ~~> around the box
            for j in range(self.ny1d):
                ipoin = j
                self.slf2d.ipob3[j] = ipoin
            for j in range(self.ny1d):
                ipoin = self.ny1d + self.nx1d - 2 + j
                self.slf2d.ipob3[self.ny1d*self.nx1d-j-1] = ipoin
            for i in range(1, self.nx1d-1):
                ipoin = i + self.ny1d - 1
                self.slf2d.ipob3[(i+1)*self.ny1d-1] = ipoin
            for i in range(1, self.nx1d-1):
                ipoin = self.nx1d + 2 * self.ny1d + i - 3
                self.slf2d.ipob3[self.ny1d*self.nx1d-i*self.ny1d-self.ny1d] = \
                    ipoin

        # ~~> Boundary points
        self.slf2d.iparam = [0, 0, 0, 0, 0, 0, 0,
                             2*self.nx1d+2*(self.ny1d-2), 0, 1]

    def put_geometry(self, file_name):

        print('   +> writing up the geometry file')

        self.slf2d.fole = {}
        self.slf2d.fole.update({'hook': open(file_name, 'wb')})
        self.slf2d.fole.update({'name': file_name})
        self.slf2d.fole.update({'endian': ">"})     # big endian
        self.slf2d.fole.update({'float': ('f', 4)})  # single precision

        self.slf2d.varnames = ['RANGE          ']
        self.slf2d.varunits = ['UI             ']
        self.slf2d.nbv1 = len(self.slf2d.varnames)
        self.slf2d.nvar = self.slf2d.nbv1
        self.slf2d.varindex = range(self.slf2d.nvar)

        print('       - Write Selafin header')
        self.slf2d.append_header_slf()

        # ~~> A few more number and a spectral template for input/output
        grbs = pygrib.open(self.dataset[0])
        for grb in grbs:
            nb_direct = grb.numberOfDirections
            nb_freq = grb.numberOfFrequencies
            break
        grbs.close()
        spec = np.zeros((nb_direct, nb_freq, self.nx1d, self.ny1d),
                        dtype=np.float)
        var = np.zeros((self.nx1d, self.ny1d), dtype=np.float)

        print('       - Write Selafin core')
        ibar = 0
        pbar = ProgressBar(maxval=len(self.slf2d.tags['times'])).start()
        for itime in range(len(self.slf2d.tags['times'])):

            self.slf2d.append_core_time_slf(self.slf2d.tags['times'][itime])
            grbs = pygrib.open(self.dataset[itime])
            for grb in grbs:
                i_i = 0
                data = grb.values.data
                data[np.where(np.absolute(data) <= 0.001)] = np.nan
                data[np.where(data == self.missing_value)] = np.nan
                data = 10. ** data
                data[np.isnan(data)] = 0.
                for i_y in range(len(self.masky)):
                    if self.masky[i_y]:
                        spec[grb.directionNumber-1,
                             grb.frequencyNumber-1, :, i_i] = \
                              data[i_y][self.maskx]
                        i_i += 1
            grbs.close()

            ibar += 1
            pbar.update(ibar)
            for i_x in range(self.nx1d):
                for i_y in range(self.ny1d):
                    var[i_x, i_y] = max(spec[:, :, i_x, i_y].ravel()) -\
                                      min(spec[:, :, i_x, i_y].ravel())
            self.slf2d.append_core_vars_slf([var.ravel()])

        pbar.finish()
        self.slf2d.fole['hook'].close()

    def set_spectral(self):

        print('   +> reseting the header of the spectral file')

        print('      - read the spectra definition')
        grbs = pygrib.open(self.dataset[0])
        for grb in grbs:
            self.nb_direct = grb.numberOfDirections
            self.nb_freq = grb.nb_freq
            self.freq = np.asarray(grb.scaledFrequencies, dtype=np.float) / \
                grb.frequencyScalingFactor
            #  /!? only so that TOMAWAC works
            self.dirc = np.asarray(grb.scaledDirections, dtype=np.float) / \
                grb.directionScalingFactor - 7.5
            break
        grbs.close()

        # ~~> sizes (spectral numbers)
        self.slf2d.nplan = 1
        self.slf2d.ndp2 = 4
        self.slf2d.ndp3 = self.slf2d.ndp2
        self.slf2d.npoin2 = self.nb_direct * self.nb_freq
        self.slf2d.npoin3 = self.slf2d.npoin2
        self.slf2d.nelem2 = self.nb_direct * (self.nb_freq-1)
        self.slf2d.nelem3 = self.slf2d.nelem2
        self.slf2d.nptfr = 2*self.nb_direct
        self.slf2d.iparam = [0, 0, 0, 0, 0, 0, 0, 2*self.nb_direct, 0, 1]

        # ~~> 2D grid (spectral grid) - TODO: use numpy here !
        self.slf2d.meshx = np.zeros(self.slf2d.npoin2, dtype=np.float)
        self.slf2d.meshy = np.zeros(self.slf2d.npoin2, dtype=np.float)
        print('      - set the mesh')
        ipoin = 0
        pbar = ProgressBar(maxval=self.slf2d.npoin2).start()
        for j_f in range(self.nb_freq):
            for i_i in range(self.nb_direct):
                self.slf2d.meshx[i_i+self.nb_direct*j_f] = \
                    self.freq[j_f]*math.sin(math.pi*self.dirc[i_i]/180.)
                self.slf2d.meshy[i_i+self.nb_direct*j_f] = \
                    self.freq[j_f]*math.cos(math.pi*self.dirc[i_i]/180.)
                ipoin += 1
                pbar.update(ipoin)
        pbar.finish()

        # ~~> Connectivity - TODO: use numpy here !
        print('      - set the connectivity')
        ielem = 0
        pbar = ProgressBar(maxval=self.slf2d.nelem3).start()
        self.slf2d.ikle3 = np.zeros((self.slf2d.nelem3, self.slf2d.ndp3),
                                    dtype=np.int)
        for j_f in range(self.nb_freq-1):
            for i_i in range(self.nb_direct):
                self.slf2d.ikle3[ielem][0] = (i_i+1) % self.nb_direct + \
                                                      j_f*self.nb_direct
                ielem += 1
        for ielem in range(self.slf2d.nelem3):
            self.slf2d.ikle3[ielem][1] = ielem
            self.slf2d.ikle3[ielem][2] = ielem + self.nb_direct
            self.slf2d.ikle3[ielem][3] = self.slf2d.ikle3[ielem][0] + \
                self.nb_direct
            pbar.update(ielem)
        pbar.finish()

        # ~~> Boundaries - TODO: use numpy here !
        pbar = ProgressBar(maxval=self.nx1d+self.ny1d).start()
        self.slf2d.ipob3 = np.zeros(self.slf2d.npoin3, dtype=np.int)
        # ~~> along the ?-axis
        for i_i in range(self.nb_direct):
            self.slf2d.ipob3[i_i] = i_i
        for i_i in range(self.nb_direct, 2*self.nb_direct):
            self.slf2d.ipob3[i_i] = self.nb_direct * \
                                          (self.nb_freq+1) - i_i
        pbar.finish()

    def append_header_grib(self):

        self.slf2d.varnames = []
        self.slf2d.varunits = []
        if self.typ == 'wave':
            # TODO: codes for waves
            raise TelemacException('... waves, not coded yet')
        elif self.typ == 'oper':
            for i in self.variables:
                if 151 == i:
                    self.slf2d.varnames.append('SURFACE PRESSURE')
                    self.slf2d.varunits.append('UI              ')
                if 165 == i:
                    self.slf2d.varnames.append('WIND VELOCITY U ')
                    self.slf2d.varunits.append('M/S             ')
                if 166 == i:
                    self.slf2d.varnames.append('WIND VELOCITY V ')
                    self.slf2d.varunits.append('M/S             ')
                if 167 == i:
                    self.slf2d.varnames.append('AIR TEMPERATURE ')
                    self.slf2d.varunits.append('DEGREES         ')
            for var in self.slf2d.varnames:
                print('    - ', var)
        elif self.typ == 'spec':
            if 251 in self.variables:
                for i in range(self.nx1d * self.ny1d):
                    self.slf2d.varnames.append(
                              ('F PT '+str(i+1)+'                ')[:16])
                    self.slf2d.varunits.append('UI              ')
            print('    - from ', self.slf2d.varnames[0], ' to ',
                  self.slf2d.varnames[-1])
        if self.slf2d.varnames == []:
            raise TelemacException(
               '... could not match requested valiable with type of file.\n\n')
        self.slf2d.nbv1 = len(self.slf2d.varnames)
        self.slf2d.nvar = self.slf2d.nbv1
        self.slf2d.varindex = range(self.slf2d.nvar)

        self.slf2d.append_header_slf()

    def append_core_time_grib(self, itime):

        self.slf2d.append_core_time_slf(self.slf2d.tags['times'][itime])

    def append_core_vars_grib(self, itime):

        if self.typ == 'wave':
            pass
            # ~~> WAVE HEIGHT == 'swh'
            # ~~> SIGNIFICANT WAVE PERIOD == 'mwp'
            # ~~> MEAN WAVE DIRECTION == 'mwd'

        elif self.typ == 'oper':
            var2d = np.zeros((self.slf2d.nvar, self.slf2d.npoin2),
                             dtype=np.float)
            grbs = pygrib.open(self.dataset[itime])
            for grb in grbs:
                if grb.indicatorOfParameter in self.variables:
                    jvar = self.variables.index(grb.indicatorOfParameter)
                    var2d[jvar, :] = np.ravel(grb.values.T)
            grbs.close()
            for jvar in range(self.slf2d.nvar):
                self.slf2d.append_core_vars_slf([var2d[jvar, :]])

        elif self.typ == 'spec':

            spec = np.zeros((self.nb_direct, self.nb_freq,
                             self.nx1d, self.ny1d),
                            dtype=np.float)
            grbs = pygrib.open(self.dataset[itime])
            ibar = 0
            maxval = self.nb_direct*self.nb_freq
            pbar = ProgressBar(maxval=maxval).start()
            for grb in grbs:
                i_i = 0
                data = grb.values.data
                data[np.where(np.absolute(data) <= 0.001)] = np.nan
                data[np.where(data == self.missing_value)] = np.nan
                data = 10. ** data
                data[np.isnan(data)] = 0.
                for i_y in range(len(self.masky)):
                    if self.masky[i_y]:
                        spec[grb.directionNumber-1, grb.frequencyNumber-1, :,
                             i_i] = \
                                  data[i_y][self.maskx]
                        i_i += 1
                ibar += 1
                pbar.update(ibar)
            pbar.finish()
            grbs.close()

            for i_x in range(self.nx1d):
                for i_y in range(self.ny1d):
                    self.slf2d.append_core_vars_slf([np.ravel(
                        spec[:, :, i_x, i_y].T)])

    def put_content(self, file_name, showbar=True):

        self.open_grib(file_name)

        print('     +> Write Selafin header')
        self.append_header_grib()

        print('     +> Write Selafin core')
        if showbar:
            pbar = ProgressBar(maxval=len(self.dataset)).start()
        for itime in range(len(self.dataset)):
            seconds = int(self.slf2d.tags['times'][itime])
            date = (datetime(*self.slf2d.datetime) +
                    timedelta(seconds=seconds)).timetuple()[0:6]
            print("        - {}-{}-{} {}:{}:{}".format(date[2], date[1],
                                                       date[0], date[3],
                                                       date[4], date[5]))
            self.append_core_time_grib(itime)
            self.append_core_vars_grib(itime)
            if showbar:
                pbar.update(itime)
        if showbar:
            pbar.finish()

        self.close_grib()
