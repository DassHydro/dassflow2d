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
from data_manip.formats.selafin import Selafin
from utils.progressbar import ProgressBar
from pretel.meshes import subdivide_mesh4

class SubSelafin(Selafin): # TODO with 3D

    def __init__(self, f):
        Selafin.__init__(self, f)
        self.ikle2, self.meshx, self.meshy, self.ipob2, \
                  self.interp, self.interp3 = \
                        subdivide_mesh4(self.ikle2, self.meshx, self.meshy)

    def put_content(self, file_name, showbar=True):
        # ~~> Doubling the number of nplan
        nplo = self.nplan
        if self.nplan > 1:
            self.nplan = 2 * self.nplan - 1
        if self.iparam[6] > 1:
            self.iparam[6] = self.nplan
        # ~~> Getting the new size of npoin2 from meshx
        np2o = self.npoin2
        self.npoin2 = len(self.meshx)
        np2n = self.npoin2
        # ~~> Setting the new size of npoin3
        np3o = self.npoin3
        self.npoin3 = self.nplan * self.npoin2
        np3n = self.npoin3
        # ~~> Getting the new size of nelem2 from self.ikle2
        self.nelem2 = len(self.ikle2)
        # ~~> Setting the new size of nelem3
        if self.nplan > 1:
            self.nelem3 = (self.nplan-1)*self.nelem2
        else:
            self.nelem3 = self.nelem2
        # ~~> Connecting
        if self.nplan > 1:
            self.ipob3 = np.ravel(np.add(np.repeat(self.ipob2, self.nplan)\
                                           .reshape((self.npoin2, self.nplan)),
                                         self.npoin2*np.arange(self.nplan)).T)
            self.ikle3 = np.repeat(self.npoin2*np.arange(self.nplan-1),
                                   self.nelem2*self.ndp3)\
                           .reshape((self.nelem2*(self.nplan-1), self.ndp3)) + \
                np.tile(np.add(np.tile(self.ikle2, 2),
                               np.repeat(self.npoin2*np.arange(2), self.ndp2)),
                        (self.nplan-1, 1))
        else:
            self.ipob3 = self.ipob2
            self.ikle3 = self.ikle2
        # ~~> Filing
        self.fole.update({'hook':open(file_name, 'wb')})
        self.fole['name'] = file_name
        self.append_header_slf()
        pbar = ProgressBar(maxval=len(self.tags['times'])).start()
        # ~~> Time stepping
        varx = np.zeros((self.nvar, self.npoin3), np.float32)
        for itime in range(len(self.tags['times'])):
            self.append_core_time_slf(itime)
            self.npoin3 = np3o           #\
            vrs = self.get_values(itime)     #|+ game of shadows
            self.npoin3 = np3n           #/
            # TODO:(JPC), convert to numpy calculations
            for ivar in range(self.nvar):
                for iplan in range(nplo):
                    nsize = 2*iplan*np2n
                    osize = iplan*np2o
                    varx[ivar][0+nsize:np2o+nsize] = \
                            vrs[ivar][0+osize:np2o+osize]
                    varx[ivar][np2o+nsize:np2o+nsize+np2n-np2o] = \
                              np.sum(vrs[ivar][self.interp+osize], axis=1)/2.
                for iplan in range(nplo-1):
                    varx[ivar][(2*iplan+1)*np2n:(2*iplan+2)*np2n] = \
                              (varx[ivar][2*iplan*np2n:(2*iplan+1)*np2n]+\
                               varx[ivar][(2*iplan+2)*np2n:(2*iplan+3)*np2n])/2.
            self.append_core_vars_slf(varx)
            pbar.update(itime)
        pbar.finish()
        self.fole['hook'].close()
