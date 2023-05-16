# -*- coding: utf-8 -*-
"""
ClassMod2D
"""

from __future__ import print_function

import copy
import json
import os
import sys
import numpy as np
from telapy.api.t2d import Telemac2d

from telapy.coupling.class_comm_2d import ClassComm2D
from telapy.coupling.class_convergence import ClassConvergence
from telapy.coupling.class_cpl_2d import ClassCpl2D


class ClassMod2D:
    """
    Relevant variables and operations to define a coupled 2D model
    """

    def __init__(self):
        """
        Constructor

        1. Initialize the coupling parallel context under
        """
        self.vars_2d = None
        self.vars_1d = None
        self.comm_2d = None
        self.conlim_2d = np.array([])
        self.conlim_co = np.array([])

        self.coupler = ClassComm2D(self)
        # Creation of the instance Telemac2d
        self.t2d = Telemac2d('T2DCAS', comm=self.coupler.comm_2d)
        self.rank = self.coupler.rank
        if self.rank == 0:
            progress = '[' + 20 * ' ' + '] ' + str(0).rjust(3) + '%'
            print(progress, file=sys.stderr, end='\r')

        self.ncsize = self.coupler.ncsize
        self.t2d.set_case()
        self.t2d.init_state_default()
        # Variable application
        self.stp = 0
        self.time = 0
        self.iter = 0
        self.converged = False
        self.inirun = False
        self.freq_res = -1
        self.maxiter = -1
        self.cplsteps = -1
        self.get_param_cpl()
        self.coupler.init_comms()
        self.dtime = self.dt_2d
        self.dtime1d = self.dt_1d
        self.nit = int(self.freq / self.dtime)
        self.ptinig_ori = 0
        self.cond1 = True
        self.l_t = 0
        self.atsave = 0
        self.ltsave = 0
        self.jsonbc_in = \
            'bc2D_restart_' + self.config + '_' + str(int(self.t_ini_2d)) \
            + '.json'
        self.jsonbc_out = \
            'bc2D_restart_' + self.config + '_' + str(int(self.t_end_2d)) \
            + '.json'
        self.fct_cpl = ClassCpl2D(self.t2d, self.coupler,
                                  self.nb_model_1d, self.nb_criteria,
                                  self.id_fr_co_2d)

        self.has_fr_co = self.fct_cpl.has_fr_co
        self.dbgcom = False
        self.last_from1d = np.zeros(self.nb_model_1d, dtype=np.float64)

        self.init_cpl()

    def get_param_cpl(self):
        """
        Get coupling parameters
        """

        try:
            with open('CurrentRunDef.json') as cfgfile:
                runcfg = json.load(cfgfile)
        except(ValueError, IOError):
            raise ValueError('Current Run Def isn\'t a correct json file.\n'
                             'Validate it with https://jsonlint.com\nSTOP')
        self.jsoncfg = runcfg['coupling_def']

        cpldef = self.jsoncfg["Coupling"]
        mods1d = self.jsoncfg["1D"]
        mod2d = self.jsoncfg["2D"][list(self.jsoncfg["2D"].keys())[0]]
        itface = self.jsoncfg["Interfaces"]
        self.inirun = runcfg["Run"]["InitialRun"]

        self.config = list(self.jsoncfg["2D"].keys())[0]
        self.t_ini_2d = runcfg["Run"]["StartTime"]
        self.t_end_2d = self.t_ini_2d + runcfg["Run"]["RunLength"]
        self.nb_model_1d = len(mods1d)
        self.cplmethod = cpldef["Method"].lower()
        if "CplStepRestart2D" in cpldef:
            self.cpl2dstart = cpldef["CplStepRestart2D"].lower()
        else:
            self.cpl2dstart = 'persist1d'
        self.maxiter = cpldef["MaxIter"]
        self.freq = cpldef["TimeStep"]
        self.dt_2d = mod2d["TimeStep"]
        self.freq_res = max(int(mod2d["OutputFreq"] / self.freq), 1)
        self.cplsteps = runcfg["Run"]["CplTSteps"]
        try:
            self.nb_site_2d = len(mod2d["OutputSites"])
        except KeyError:
            self.nb_site_2d = 0
        if self.nb_site_2d > 0.:
            self.tab_site_2d = np.array(mod2d["OutputSites"])
        self.crit_arret_h = [i["ConvCriteria"]["Height"] for i in itface]
        self.crit_arret_v = [i["ConvCriteria"]["Velocity"] for i in itface]
        self.nit = int(self.freq / self.dt_2d)
        self.nb_criteria = 6
        # WARNING: we rely for the moment on the hypothesys
        # of common 1d timestepping
        self.dt_1d = mods1d[list(mods1d.keys())[0]]["TimeStep"]
        self.nit_1d = int(self.freq / self.dt_1d) + 1
        self.type_cl_in_1d = \
            np.array([1 if i["Condition1D"].lower() ==
                     "waterlevel" else 2 for i in itface])
        self.pos_model_1d = \
            np.array([1 if i["1DPosition"].lower() ==
                     "upstream" else 2 for i in itface])
        self.id_fr_co_2d = np.array([i["LiqBdry2D"] for i in itface])

    @staticmethod
    def checkfile(cas_file):
        """
        Check if file exists
        @param cas_file (str):  path file
        """
        if not os.path.exists(cas_file):
            return False
        return True

    def init_cpl(self):
        """
         Initialization of coupling model variables
         @return nothing
        """

        self.t2d.set('MODEL.DT', self.dtime)
        self.t2d.set('MODEL.NTIMESTEPS', self.nit)

        self.t2d.save_state()
        cond_json = False
        cote = None
        deb = None
        if self.cplmethod == 'additiveschwarz' and self.has_fr_co:
            cote = self.t2d.get_array('MODEL.COTE')
            deb = self.t2d.get_array('MODEL.DEBIT')

            if self.checkfile(self.jsonbc_in):
                cond_json = True
                with open(self.jsonbc_in) as json_data:
                    dico = json.load(json_data)
                for id_model in dico['model']:
                    idm = int(id_model)
                    if self.type_cl_in_1d[idm] == 1:
                        cote[self.id_fr_co_2d[idm] - 1] = \
                            float(dico['conlim_co'][idm])
                    elif self.type_cl_in_1d[idm] == 2:
                        deb[self.id_fr_co_2d[idm] - 1] = \
                            float(dico['conlim_co'][idm])
                self.last_from1d = dico['last_from1d']

                self.t2d.set_array('MODEL.COTE', cote)
                self.t2d.set_array('MODEL.DEBIT', deb)
                # os.remove(self.jsonbc_in)

        # AP insert something here to store the values
        # at the coupling interfaces
        # AP as read in from the model definition
        # (unless self.cpl2dstart 'persist2d')
        if self.cplmethod == 'additiveschwarz' \
                and self.cpl2dstart != 'persist2d' \
                and self.has_fr_co:

            self.ini_frliq = np.zeros(self.nb_model_1d, dtype=np.float64)

            for i, typc in enumerate(self.type_cl_in_1d):
                # ! cas 1 :  The 1D model have Q in BC, it send  Z
                if typc == 1:
                    self.ini_frliq[i] = \
                        copy.deepcopy(cote[self.id_fr_co_2d[i] - 1])

                # ! cas 2 :  The 1D model have Z in BC, it send  Q
                elif typc == 2:
                    self.ini_frliq[i] = \
                        copy.deepcopy(deb[self.id_fr_co_2d[i] - 1])
            if not cond_json:
                self.last_from1d = copy.deepcopy(self.ini_frliq)

            if self.dbgcom and self.rank == 0:
                print('DBG add pr1d or pe1d init_step ini_frliq\n',
                      self.ini_frliq)

        self.conv = ClassConvergence(self.nb_model_1d, self.maxiter,
                                     output=self.rank == 0)
        self.conv.set_criteria(height=self.crit_arret_h,
                               velocity=self.crit_arret_v)

        if self.rank == 0:
            print("FIN INITIALISATION DES PARAMETRE COUPLAGE")

    def run(self):
        """
        Integrate T2D on a full coupling step
        @return nothing
        """
        if self.rank == 0:
            print('DBG run: step', self.stp, 'iter', self.iter)
        self.init_step()
        while self.l_t < (self.stp * self.nit + self.nit):
            self.step()
        self.time = (self.stp + 1) * self.freq

    def init_step(self):
        """
        Save the model state, get the boundary conditions,
        and treatment of the boundary conditions
        @return nothing
        """

        if self.iter == 0:
            self.t2d.save_state()
            self.atsave = self.t2d.get('MODEL.AT')
            self.ltsave = self.t2d.get('MODEL.LT')
        else:
            self.t2d.restore_state()
            self.t2d.set('MODEL.AT', self.atsave)
            self.t2d.set('MODEL.LT', self.ltsave)

        # AP insert here the handling of the initial bc update
        # accordingly to method
        # AP and (id additive) restart variants self.cpl2dstart
        if self.cplmethod == 'multiplicativeschwarz' and self.has_fr_co:
            self.coupler.get_cl()
            self.fct_cpl.interp_conlim(self.nit, self.dt_2d,
                                       self.dt_1d, self.conlim_co)
        elif self.cplmethod == 'additiveschwarz' and self.has_fr_co:
            if self.cpl2dstart == 'previous1d':
                self.update_previous1d()
            elif self.cpl2dstart == 'persist1d':
                self.update_persist1d()

            elif self.cpl2dstart == 'persist2d':
                self.update_persist2d()
        self.l_t = self.t2d.get('MODEL.LT')
        self.conlim_2d[:, self.l_t - self.stp * self.nit] = \
            self.fct_cpl.comput_conlim(self.pos_model_1d,
                                       self.type_cl_in_1d)

    def update_previous1d(self):
        """
        Update boundary conditions with previous1d method:
        start with the last conditions received from the 1D model on
        the previous coupling step and interpolate them in time.
        """
        if self.stp == 0 and self.iter == 0:
            self.fct_cpl.interp_conlim_co = \
                np.repeat([self.ini_frliq], self.nit + 1, axis=0)
        else:
            self.coupler.get_cl()
            self.fct_cpl.interp_conlim(self.nit, self.dt_2d,
                                       self.dt_1d, self.conlim_co)
        if self.dbgcom and self.rank == 0:
            print('DBG add pr1d conlim_co\n', self.fct_cpl.interp_conlim_co)

    def update_persist1d(self):
        """
        Update boundary conditions with persist1d method:
        start with the last 1d state values at the coupling
        sections (constant in time)
        """
        if self.stp == 0 and self.iter == 0:
            if self.inirun:
                self.fct_cpl.interp_conlim_co = \
                    np.repeat([self.ini_frliq], self.nit + 1, axis=0)
                self.last_from1d = copy.deepcopy(self.ini_frliq)
            else:
                self.last_from1d = (self.ini_frliq + self.last_from1d) * 0.5
                self.fct_cpl.interp_conlim_co = \
                    np.repeat([self.last_from1d], self.nit + 1, axis=0)

            if self.dbgcom and self.rank == 0:
                print('DBG add pe1d s0 i0 conlim_co\n',
                      self.fct_cpl.interp_conlim_co)
                print('DBG add pe1d s0 i0 last_from1d', self.last_from1d)
        elif self.stp > 0 and self.iter == 0:
            self.coupler.get_cl()
            if self.dbgcom and self.rank == 0:
                print('DBG add pe1d s>0 i0 ori last_from1d', self.last_from1d)
                print('DBG add pe1d s>0 i0 received co', self.conlim_co[:, -1])
            self.last_from1d = (self.conlim_co[:, -1] + self.last_from1d) * 0.5
            self.fct_cpl.interp_conlim_co = \
                np.repeat([self.last_from1d], self.nit + 1, axis=0)
            if self.dbgcom and self.rank == 0:
                print('DBG add pe1d s>0 i0 new last_from1d', self.last_from1d)
                print('DBG add pe1d s>0 i0 conlim_co\n',
                      self.fct_cpl.interp_conlim_co)
        else:
            self.coupler.get_cl()
            self.last_from1d = copy.deepcopy(self.conlim_co[:, -1])
            self.fct_cpl.interp_conlim(self.nit, self.dt_2d,
                                       self.dt_1d, self.conlim_co)
            if self.dbgcom and self.rank == 0:
                print('DBG add pe1d s>0 || i>0 last_from1d', self.last_from1d)
                print('DBG add pe1d s>0 || i>0 conlim_co\n',
                      self.fct_cpl.interp_conlim_co)

    def update_persist2d(self):
        """
        Update boundary conditions with persist2d method:
        Start with the last conditions received from the 2D model
            on the previous
        coupling step but keep only the last time step value (constant in time)
        """

        if self.iter == 0:
            if self.stp > 0:
                self.coupler.get_cl()
            self.fct_cpl.interp_conlim_co = \
                np.repeat([self.ini_frliq], self.nit + 1, axis=0)
        else:
            self.coupler.get_cl()
            self.fct_cpl.interp_conlim(self.nit, self.dt_2d,
                                       self.dt_1d, self.conlim_co)
        if self.dbgcom and self.rank == 0:
            print('DBG add pe2d conlim_co\n', self.fct_cpl.interp_conlim_co)

    def step(self):
        """
        Change boundary conditions and compute time step
        @return nothing
        """

        if self.has_fr_co:
            self.fct_cpl.cl_modif(self.type_cl_in_1d,
                                  self.fct_cpl.interp_conlim_co[(self.l_t + 1 -
                                                                 self.stp *
                                                                 self.nit), :])

        self.t2d.run_one_time_step_compute()
        self.l_t = self.t2d.get('MODEL.LT')
        self.conlim_2d[:, self.l_t - self.stp * self.nit] = \
            self.fct_cpl.comput_conlim(self.pos_model_1d,
                                       self.type_cl_in_1d)

    def transmit_step(self):
        """
        Compute the convergence variables,
        transmit the boundary conditions of 1D model,
        and get the convergence variable
        @return nothing
        """

        self.fct_cpl.var_interface(self.pos_model_1d, self.type_cl_in_1d,
                                   self.vars_2d)
        if self.rank == 0:
            self.coupler.transmit_cl_1d()
            self.coupler.get_vars_1d()
            self.converged = self.conv.main(self.iter, self.vars_1d,
                                            self.vars_2d, ib_tps=self.stp)
        if self.ncsize > 1:
            self.converged = self.coupler.comm_2d.bcast(self.converged, root=0)

    def write_output(self):
        """
        Write output
        @return nothing
        """

        step = self.stp + 1

        if self.rank == 0:
            p_c = int(step * 100. / self.cplsteps)
            pci = int((p_c - 1) / 5) + 1
            progress = '[' + pci * '=' + (20 - pci) * ' ' + '] ' + \
                str(p_c).rjust(3) + '%'
            print(progress, file=sys.stderr, end='\r')

        if step % self.freq_res == 0 and step != self.cplsteps:
            if self.rank == 0:
                print("*********** WRITE OUTPUT FOR {} S ***********"
                      .format(self.t2d.get('MODEL.AT')))
            self.t2d.run_one_time_step_res()

    def finalize(self):
        """
        Finalize the instances of Telemac 2D
        @return nothing
        """
        if self.rank == 0:
            print('\n', file=sys.stderr)
            print("*********** WRITE FINAL OUTPUT ***********")
        self.t2d.run_one_time_step_res()
        # write val cl_pc pour reprise add
        self.write_cl_restart()
        self.conv.finalize()
        self.t2d.finalize()
        # Instance delete
        del self.t2d
        if self.rank == 0:
            print('My work is done.')

    def write_cl_restart(self):
        """
        Dump last bc got from 1D
        so that restarting can be almost perfectly smooth
        @return nothing
        """
        if self.cplmethod == 'additiveschwarz' and self.has_fr_co:
            self.coupler.get_cl()

            if self.rank == 0:
                dico = {'model': [], 'conlim_co': [], 'last_from1d': []}
                for i, val in enumerate(self.conlim_co[:, -1]):
                    dico['model'].append(i)
                    dico['conlim_co'].append(val)
                    dico['last_from1d'].append(self.last_from1d[i])

                # print(dico)
                with open(self.jsonbc_out, 'w') as outfile:
                    json.dump(dico, outfile)
