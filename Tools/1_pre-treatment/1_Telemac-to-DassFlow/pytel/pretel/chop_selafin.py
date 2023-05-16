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
from data_manip.formats.parafins import Parafins
from data_manip.extraction.parser_selafin import subset_variables_slf

class ChopSelafin(Parafins):

    chop_from = 0
    chop_step = 1
    chop_stop = -1

    def __init__(self, f, times=None, vrs=None, root=None):
        Parafins.__init__(self, f, root)
        if vrs != None:
            self.update_vars(vrs)
        if times != None:
            tfrom, tstep, tstop = times
            self.update_times(tfrom, tstep, tstop)

    def update_vars(self, vrs):
        allvars = self.slf.varnames
        allvars.extend(self.slf.cldnames)
        allunits = self.slf.varunits
        allunits.extend(self.slf.cldunits)
        self.slf.varindex = subset_variables_slf(vrs, allvars)[0]
        nbv1 = 0
        varnames = []
        varunits = []
        nbv2 = 0
        cldnames = []
        cldunits = []
        for i in self.slf.varindex:
            if i <= self.slf.nbv1:
                nbv1 += 1
                varnames.append(allvars[i])
                varunits.append(allunits[i])
            else:
                nbv2 += 1
                varnames.append(allvars[i])
                varunits.append(allunits[i])
        self.slf.nbv1 = nbv1
        self.slf.varnames = varnames
        self.slf.varunits = varunits
        self.slf.nbv2 = nbv2
        self.slf.cldnames = cldnames
        self.slf.cldunits = cldunits

    def update_times(self, tfrom, tstep, tstop): # /!\ starts from 0!
        if tfrom > 0:
            tfrom -= 1
        if tstop > 0:
            tstop -= 1
        ltime = len(self.slf.tags['times'])
        if tstop < 0:
            tstop = max(0, ltime+tstop)
        if tfrom < 0:
            tfrom = max(0, ltime+tfrom)
        tfrom = min(ltime-1, tfrom)
        tstop = max(min(ltime-1, tstop), tfrom)
        self.chop_from = tfrom
        self.chop_step = max(1, tstep)
        self.chop_stop = tstop
        att = []
        ats = []
        for itime in \
                range(ltime)[self.chop_from:self.chop_stop+1:self.chop_step]:
            ats.append(self.slf.tags['times'][itime])
            att.append(self.slf.tags['cores'][itime])
        self.slf.tags['times'] = np.asarray(ats)
        self.slf.tags['cores'] = att

