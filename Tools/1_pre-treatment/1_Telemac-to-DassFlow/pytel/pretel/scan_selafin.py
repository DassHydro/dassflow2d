r"""@author Sebastien E. Bourban

"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import numpy as np
# ~~> dependencies towards other modules
# ~~> dependencies towards other modules
from pretel.chop_selafin import ChopSelafin

# /!\ does not support Parafins yet -- because of the print order of print Core
class ScanSelafin(ChopSelafin):

    def __init__(self, f, times=None, vrs=None):
        ChopSelafin.__init__(self, f)
        if vrs != None:
            ChopSelafin.update_vars(self, vrs)
        if times != None:
            tfrom, tstep, tstop = times
            ChopSelafin.update_times(self, tfrom, tstep, tstop)

    def print_header(self):
        if self.slf.file['endian'] == ">":
            print('This file appears to be coded in "big endian"')
        else:
            print('This Selafin file appears to be coded in "little endian"')
        if self.slf.file['float'] == ('f', 4):
            print(5*' '+'and the floats are assumed to be SINGLE PRECISION\n')
        else:
            print(5*' '+'and the floats are assumed to be DOUBLE PRECISION\n')
        print("title        :    <{}>".format(self.slf.title))
        if self.slf.iparam[9] == 1:    # /!\ needs proper formatting
            print("Date / Time  : {}-{}-{} {}:{}:{}"\
                .format(self.slf.datetime[2],
                        self.slf.datetime[1],
                        self.slf.datetime[0],
                        self.slf.datetime[3],
                        self.slf.datetime[4],
                        self.slf.datetime[5]))
        else:
            print("No date in file")
        if len(self.slf.varnames) > 0:
            print("Variables    :\n   - " + \
                    "\n   - ".join(v+u for v, u in zip(self.slf.varnames,
                                                       self.slf.varunits)))
        if len(self.slf.cldnames) > 0:
            print("Clandestines :\n   - " +
                  "\n   - ".join(v+u for v, u in zip(self.slf.cldnames,
                                                     self.slf.cldunits)))
        print("NUMBERs      :")
        print("   - nplan* = {}".format(self.slf.iparam[6]))
        print("   - nptfr* = {}".format(self.slf.iparam[7]))
        print("   - iface* = {}".format(self.slf.iparam[8]))
        print("   - nelem3 = {}".format(self.slf.nelem3))
        print("   - npoin3 = {}".format(self.slf.npoin3))
        print("   - ndp3   = {}".format(self.slf.ndp3))
        print("   - nplan  = {}".format(self.slf.nplan))
        print("   - x_orig/y_orig = {}, {}".format(self.slf.iparam[2],
                                              self.slf.iparam[3]))
        if self.slf.nplan > 1:
            print("   - nelem2 = {}".format(self.slf.nelem2))
            print("   - npoin2 = {}".format(self.slf.npoin2))
            print("   - ndp2   = {}".format(self.slf.ndp2))
        if self.slf.npoin2 > 0:
            x_min = np.min(self.slf.meshx)
            y_min = np.min(self.slf.meshy)
            x_max = np.max(self.slf.meshx)
            y_max = np.max(self.slf.meshy)
            print("MESH         : / min: [ {};{} ] / max: [ {};{} ]".format(\
                    x_min, y_min, x_max, y_max))
        print("ARRAYs       :")
        print("   - ikle  : / min: [ {} ]  / max: [ {} ] {}"\
              .format(np.min(self.slf.ikle3)+1, np.max(self.slf.ikle3)+1,
                      self.slf.ikle3+1))
        print("   - ipobo  : / min: [ {} ]  / max: [ {} ] {}"\
              .format(np.min(self.slf.ipob3), np.max(self.slf.ipob3),
                      self.slf.ipob3))

    def print_core(self):
        for var in range(self.slf.nbv1):
            print("Variable     : "+self.slf.varnames[var])
            for itime in range(len(self.slf.tags['times'])):
                varsor = self.slf.get_variables_at(itime,
                                                   [self.slf.varindex[var]])
                if self.slf.npoin2 > 0:
                    print("    / Time: {} / min: {} / max: {}".format(\
                          self.slf.tags['times'][itime],
                          np.min(varsor[0]), np.max(varsor[0])))
        for var in range(self.slf.nbv2):
            print("Clandestine  : "+self.slf.cldnames[var])
            for itime in range(len(self.slf.tags['times'])):
                idx = var+self.slf.nbv1
                varsor = self.slf.get_variables_at(itime,
                                                   [self.slf.varindex[idx]])
                if self.slf.npoin2 > 0:
                    print("    / Time: {} / min: {} / max: {}".format(\
                          self.slf.tags['times'][itime],
                          np.min(varsor[0]), np.max(varsor[0])))

    def print_time_summary(self):
        print("Number of times : {}".format(len(self.slf.tags['times'])))
        print("First time step : {}".format(self.slf.tags['times'][0]))
        print("Last time step  : {}".format(\
                self.slf.tags['times'][-1]))
        if len(self.slf.tags['times']) > 1:
            times = self.slf.tags['times']
            print("Time step       : "+ \
            str((times[-1] - times[0]) / (len(times)-1)))
        else:
            print("Only one time frame")

