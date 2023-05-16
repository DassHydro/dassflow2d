r"""@author Sebastien E. Bourban

"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards other modules
from pretel.alter_selafin import AlterSelafin
from utils.progressbar import ProgressBar
# ~~> dependencies towards standard python
import numpy as np

class TransfSelafin(AlterSelafin):

    def __init__(self, f, times=None, vrs=None, root=None, points=None):
        AlterSelafin.__init__(self, f, root)
        self.calcs = []
        # ~~> history has no impact on potential resampling of
        #     self.slf.tags['times']
        if times == (1, 1, -1):
            self.history = np.arange(2006.07, 2108.75, 0.34, dtype=float)
        else:
            tfrom, tstep, tstop = times
            self.history = np.arange(tfrom, tstop, tstep, dtype=float)
        # ~~> extraction restricted to points
        if points != None:
            pass

    def calc_free_surface_from_artemis(self):
        # ~~> Dependencies
        #     needs all the components to recreate the free surface signal
        args = range(self.slf.nbv1)
        # ~~> New variable name
        calcs = {'vars':[["WAVE SURFACE    ", "M               "]]}
        # ~~> Initial value 0
        def init(vrs, ivars, t_0, t_i):
            return [np.zeros(self.slf.npoin3, dtype=np.float64)]
        calcs.update({'init':(init, args)})
        # ~~> Computation of free surface
        def calc(vrs, ivars, ttime, t_i, lastt, vari):
            # work in double precision
            vars64 = np.asarray(vrs, dtype=np.float64)
            sqrt2 = np.sqrt(2.0)
            for dir_n in range(len(vars64)//2):
                # dfr is an artefact to avoid phase locking
                dfr = np.float64(0.04/t_i)*\
                        np.float64(dir_n+1-(len(vars64)//2+1)/2)
                ivars1 = ivars[2*dir_n]
                ivars2 = ivars[2*dir_n+1]
                coeff = -2.0*np.pi*np.float64(1.0/t_i+dfr)*np.float64(ttime)
                vari[0] += vars64[ivars1]*np.cos(coeff+vars64[ivars2])\
                           /(2.0*sqrt2)
            return vari

        calcs.update({'calc':(calc, args)})
        # ~~> Conclusion step
        def stop(t_0, t_i, vari):
            return vari
        calcs.update({'stop':stop})
        # ~~> Store
        self.calcs.append(calcs)

    def put_content(self, file_name, showbar=True):
        # ~~> Output file header
        self.slf.fole.update({'hook': open(file_name, 'wb')})
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

        # ~~> Time stepping
        print('\n      > Input signal based on:')
        print('      - '+str(self.slf.nvar/2)+' direction(s)')
        print('      - for the following periods:'+repr(self.slf.tags['times']))
        # ~~> Time stepping
        print('\n      > Going through time ('+str(len(self.history))+\
                'time steps) :')
        pbar = ProgressBar(maxval=len(self.history)).start()

        for t_i in range(len(self.history)):
            ttime = self.history[t_i]
            # ~~> Initialise vari to 0 for each time step
            vari = []
            t_0 = self.slf.tags['times'][0]
            for calc, icalc in zip(self.calcs, range(len(self.calcs))):
                fct, args = calc['init']
                vari.append(fct(vari, args, t_0, t_0))
            # ~~> Sweeps through wave periods, adding up all components of
            # free surface into "vari" for a given time
            #     in this particular case self.slf.tags['times'] holds
            #     wave periods, not times
            for t_p in range(len(self.slf.tags['times'])):
                vrs = self.get_palues(t_p)
                for icalc, calc in enumerate(self.calcs):
                    fct, args = calc['calc']
                    times = self.slf.tags['times']
                    vari[icalc] = fct(vrs, args, ttime, times[t_p],
                                      times[len(times)-1], vari[icalc])
            # ~~> Print time record
            self.slf.append_core_time_slf(float(ttime))
            for icalc, calc in enumerate(self.calcs):
                fct = calc['stop']
                self.slf.append_core_vars_slf(fct(t_0, ttime, vari[icalc]))
            pbar.update(t_i)

        pbar.finish()
        self.slf.fole['hook'].close()
