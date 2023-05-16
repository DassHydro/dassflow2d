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
from utils.progressbar import ProgressBar
from data_manip.extraction.parser_selafin import subset_variables_slf

class CrunchSelafin(AlterSelafin):

    def __init__(self, f, times=None, vrs=None, root=None):
        AlterSelafin.__init__(self, f, root)
        self.calcs = []
        if times != None:
            tfrom, tstep, tstop = times
            AlterSelafin.update_times(self, tfrom, tstep, tstop)

    def calc_peak_time_modulo_m2(self):
        # ~~> Dependancies
        args = subset_variables_slf("FREE SURFACE", self.slf.varnames)[0]
        # ~~> New variable name
        calcs = {'vars':[["TIME OF PEAK    ", "H               "],
                         ["SURFACE AT PEAK ", "M               "]]}
        # ~~> Initial value for new variable
        def init(vrs, ivars, t_0, t_i):
            return [np.zeros(self.slf.npoin3, dtype=np.float32),
                    np.zeros(self.slf.npoin3, dtype=np.float32)]
        calcs.update({'init':(init, args)})
        # ~~> Calculation for new variable
        def calc(vrs, ivars, t_0, t_i, vari):
            for ipoin in range(self.slf.npoin3):
                if vrs[ivars[0]][ipoin] > vari[1][ipoin]:
                    # Modulo M2 (this could be an argument)
                    vari[0][ipoin] = ((t_i-t_0)/3600.0)%12.42
                    vari[1][ipoin] = vrs[ivars[0]][ipoin]
            return vari
        calcs.update({'calc':(calc, args)})
        # ~~> Conclusion step for new variable
        def stop(t_0, t_i, vari):
            return vari
        calcs.update({'stop':stop})
        # ~~> Store
        self.calcs.append(calcs)

    def calc_surface_range(self):
        # ~~> Dependancies
        args = subset_variables_slf("FREE SURFACE", self.slf.varnames)[0]
        # ~~> New variable name
        calcs = {'vars':[["SURFACE RANGE   ", "M               "]]}
        # ~~> Initial value for new variable
        def init(vrs, ivars, t_0, t_i):
            return [np.array(vrs[ivars[0]], copy=True),
                    np.array(vrs[ivars[0]], copy=True)]
        calcs.update({'init':(init, args)})
        # ~~> Calculation for new variable
        def calc(vrs, ivars, t_0, t_i, vari):
            return [np.minimum(vrs[ivars[0]], vari[0]),
                    np.maximum(vrs[ivars[0]], vari[1])]
        calcs.update({'calc':(calc, args)})
        # ~~> Conclusion step for new variable
        def stop(t_0, t_i, vari):
            return [vari[1]-vari[0]]
        calcs.update({'stop':stop})
        # ~~> Store
        self.calcs.append(calcs)

    def calc_residual_velocity(self):
        # ~~> Dependancies
        args = subset_variables_slf("VELOCITY U;VELOCITY V",
                                    self.slf.varnames)[0]
        # ~~> New variable name
        calcs = {'vars':[["RESIDUAL U      ", "M/S             "],
                         ["RESIDUAL V      ", "M/S             "]]}
        # ~~> Initial value for new variable
        def init(vrs, ivars, t_0, t_i):
            return [np.zeros(self.slf.npoin3, dtype=np.float32),
                    np.zeros(self.slf.npoin3, dtype=np.float32)]
        calcs.update({'init':(init, args)})
        # ~~> Calculation for new variable
        def calc(vrs, ivars, t_0, t_i, vari):
            vari[0] += vars[ivars[0]]
            vari[1] += vars[ivars[1]]
            return vari
        calcs.update({'calc':(calc, args)})
        # ~~> Conclusion step for new variable
        def stop(t_0, t_i, vari):
            return vari/(t_i-t_0)
        calcs.update({'stop':stop})
        # ~~> Store
        self.calcs.append(calcs)

    def calc_maximum_speed(self):
        # ~~> Dependancies
        args = subset_variables_slf("VELOCITY U;VELOCITY V",
                                    self.slf.varnames)[0]
        # ~~> New variable name
        calcs = {'vars':[["MAXIMUM SPEED   ", "M/S             "]]}
        # ~~> Initial value for new variable
        def init(vrs, ivars, t_0, t_i):
            return [np.zeros(self.slf.npoin3, dtype=np.float32)]
        calcs.update({'init':(init, args)})
        # ~~> Calculation for new variable
        def calc(vrs, ivars, t_0, t_i, vari):
            tmp = np.power(vrs[ivars[0]], 2)+np.power(vrs[ivars[1]], 2)
            return [np.maximum(vari[0], np.power(tmp, (1.0/2.0)))]
        calcs.update({'calc':(calc, args)})
        # ~~> Conclusion step for new variable
        def stop(t_0, t_i, vari):
            return vari
        calcs.update({'stop':stop})
        # ~~> Store
        self.calcs.append(calcs)

    def put_content(self, file_name, showbar=True):
        # ~~> Sweep through time steps, saving "vari"
        vari = []
        initiate = True
        if showbar:
            pbar = ProgressBar(maxval=len(self.slf.tags['times'])).start()
        t_0 = self.slf.tags['times'][0]
        for itime in range(len(self.slf.tags['times'])):
            vrs = self.get_palues(itime)
            if initiate:
                for icalc, calc in enumerate(self.calcs):
                    fct, args = calc['init']
                    vari.append(fct(vrs, args, t_0,
                                    self.slf.tags['times'][itime]))
                initiate = False
            else:
                for icalc, calc in enumerate(self.calcs):
                    fct, args = calc['calc']
                    vari[icalc] = fct(vrs, args, t_0,
                                      self.slf.tags['times'][itime],
                                      vari[icalc])
            if showbar:
                pbar.update(itime)
        if showbar:
            pbar.finish()
        # ~~> Header
        self.slf.fole.update({'hook':open(file_name, 'wb')})
        self.slf.fole['name'] = file_name
        self.slf.nbv1 = 0
        self.slf.nbv2 = 0
        self.slf.varnames = []
        self.slf.varunits = []
        for calc in self.calcs:
            for cname, cunit in calc['vars']:
                self.slf.varnames.append(cname)
                self.slf.varunits.append(cunit)
                self.slf.nbv1 += 1
        self.slf.append_header_slf()
        # ~~> Core
        # TODO: writing only for first timestep
        self.slf.append_core_time_slf(0)
        for icalc, calc in enumerate(self.calcs):
            fct = calc['stop']
            itime = t_0
            self.slf.append_core_vars_slf(fct(t_0, itime, vari[icalc]))
        self.slf.fole['hook'].close()
