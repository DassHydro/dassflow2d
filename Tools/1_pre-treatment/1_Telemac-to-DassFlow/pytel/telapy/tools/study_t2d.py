"""
Strucuture for a Telemac study
"""
import numpy as np
from telapy.api.t2d import Telemac2d
from mpi4py import MPI
import utils.polygon as polygon


class StudyTelemac2D(object):
    def __init__(self, studyFiles, tps_obs, poly_obs=None):
        # Initialization of the observation part
        self.tps_obs = tps_obs
        self.studyFiles = studyFiles
        self.poly_obs = poly_obs

    def h_x(self, k):
        # k : Strickler coefficient
        h, _, _ = self.run_telemac2d(k, finalize=True)
        y = np.asmatrix(np.ravel(h)).T
        return y

    def run_telemac2d(self, x, finalize=False):
        # creation of the telemac instance with compilation of the user fortran
        self.t2d = Telemac2d(self.studyFiles['t2d.cas'],
                             user_fortran=self.studyFiles['t2d.f'], stdout=0,
                             comm=MPI.COMM_WORLD)
        # Telemac Steering file reading
        self.t2d.set_case()
        # This is the initialization part
        self.t2d.init_state_default()
        self.final_time = self.t2d.get("MODEL.NTIMESTEPS")
        self.npoin = self.t2d.get('MODEL.NPOIN')
        for i in range(self.npoin):
            self.t2d.set('MODEL.CHESTR', x, i)
        self.point_obs = []
        if self.poly_obs is not None:
            for i in range(self.npoin):
                if polygon.is_in_polygon(self.t2d.get('MODEL.X', i),
                                         self.t2d.get('MODEL.Y', i),
                                         self.poly_obs):
                    self.point_obs.append(i)

        h = []
        u = []
        v = []
        if 0 in self.tps_obs:
            if len(self.point_obs) > 0:
                h.append(np.zeros((1, len(self.point_obs))))
                u.append(np.zeros((1, len(self.point_obs))))
                v.append(np.zeros((1, len(self.point_obs))))
                for k in self.point_obs:
                    h[self.tps_obs.index(0)][self.point_obs.index(k)] = \
                        self.t2d.get('MODEL.WATERDEPTH', k)
                    u[self.tps_obs.index(0)][self.point_obs.index(k)] = \
                        self.t2d.get('MODEL.VELOCITYU', k)
                    v[self.tps_obs.index(0)][self.point_obs.index(k)] = \
                        self.t2d.get('MODEL.VELOCITYV', k)
        # Beginning of the computation
        ntimesteps = self.t2d.get("MODEL.NTIMESTEPS")

        for i in range(ntimesteps):
            self.t2d.run_one_time_step()
            if i+1 in self.tps_obs:
                if len(self.point_obs) > 0:
                    h.append(np.zeros((1, len(self.point_obs))))
                    u.append(np.zeros((1, len(self.point_obs))))
                    v.append(np.zeros((1, len(self.point_obs))))
                    for k in self.point_obs:
                        h[self.tps_obs.index(i+1)][self.point_obs.index(k)] = \
                            self.t2d.get('MODEL.WATERDEPTH', k)
                        u[self.tps_obs.index(i+1)][self.point_obs.index(k)] = \
                            self.t2d.get('MODEL.VELOCITYU', k)
                        v[self.tps_obs.index(i+1)][self.point_obs.index(k)] = \
                            self.t2d.get('MODEL.VELOCITYV', k)
        if finalize is True:
            self.t2d.finalize()
            del self.t2d
        return h, u, v
