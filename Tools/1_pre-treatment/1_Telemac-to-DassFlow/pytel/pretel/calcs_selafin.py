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
from pretel.alter_selafin import AlterSelafin
from data_manip.extraction.parser_selafin import subset_variables_slf
from utils.progressbar import ProgressBar

class CalcsSelafin(AlterSelafin):

    def __init__(self, f, times=None, vrs=None, root=None):
        AlterSelafin.__init__(self, f, root)
        self.calcs = []
        if times != None:
            tfrom, tstep, tstop = times
            AlterSelafin.update_times(self, tfrom, tstep, tstop)

    def calc_water_depth(self):
        self.slf.varnames.append("WATER DEPTH     ")
        self.slf.varunits.append("M               ")
        self.slf.nbv1 += 1
        args = subset_variables_slf("FREE SURFACE;BOTTOM", self.slf.varnames)[0]
        def calc(vrs, ivars):
            return [vrs[ivars[0]]-vrs[ivars[1]]]
        self.calcs.append([calc, args])

    def calc_kinetic_energy(self):
        self.slf.varnames.append("KINETIC ENERGY  ")
        self.slf.varunits.append("?               ")
        self.slf.nbv1 += 1
        args = subset_variables_slf("VELOCITY U;VELOCITY V",
                                    self.slf.varnames)[0]
        def calc(vrs, ivars):
            var1 = np.square(vrs[ivars[0]])
            var2 = np.square(vrs[ivars[1]])
            return [np.power(var1+var2, (3.0/2.0))]
        self.calcs.append([calc, args])

    def put_content(self, file_name):
        self.slf.fole.update({'hook': open(file_name, 'wb')})
        self.slf.file['name'] = file_name
        pbar = ProgressBar(maxval=len(self.slf.tags['times'])).start()
        self.slf.append_header_slf()
        # ~~> Time stepping
        for itime in range(len(self.slf.tags['times'])):
            self.slf.append_core_time_slf(itime)
            vrs = self.get_palues(itime)
            self.slf.append_core_vars_slf(vrs)
            for fct, args in self.calcs:
                self.slf.append_core_vars_slf(fct(vrs, args))
            pbar.update(itime)
        pbar.finish()
        self.slf.fole['hook'].close()
