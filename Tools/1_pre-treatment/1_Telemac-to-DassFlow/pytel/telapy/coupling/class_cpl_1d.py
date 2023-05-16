# coding: utf-8
"""
ClassCpl1D class
================
"""
import copy
import json
import os
from math import sqrt

import numpy as np
from telapy.coupling.class_convergence import ClassConvergence


class ClassCpl1D:
    """
    Coupling tools for the 1D models
    """

    def __init__(self, model):
        """
        Constructor
            1. Receive and store the pointer to the model instance
            2. Recover the initial definition of the model timestepping
            3. Store a permanent copy of the initial timestepping
            4. Compute the index of the coupling section
        @param model (object) ClassMod1D model, pointer to 1D model instance
        """
        model.study.masc.set('Model.InitTime', model.ti_ini)
        model.tstinstp = int(model.cplfreq / model.d_t)
        model.tbcinstp = model.tstinstp + 1

        model.nbbc = model.study.masc.get_nb_cl()
        indices_names = \
            [model.study.masc.get_name_cl(k + 1) for k in range(model.nbbc)]
        idxbc = [x[0] for x in indices_names]
        namebc = [x[1] for x in indices_names]

        model.bc2 = np.zeros((model.tbcinstp, model.nbbc), dtype=np.float64)
        model.hydlaw = []

        for law in range(model.nbbc):
            hydlaw = HydLaw(model, law, namebc[law], idxbc[law] - 1)
            model.hydlaw.append(hydlaw)
            if idxbc[law] == model.coupling_law:
                model.cpllaw = model.hydlaw[law]

        if model.cpl1dstart != 'persist1d':
            # Recover the initial state variables to be used
            # as initial forcings
            if model.cpllaw.type == 1:
                model.cploutnam = 'State.Z'
                model.ini_frliq = \
                    model.study.masc.get('State.Q', model.cplsect)
            elif model.cpllaw.type == 2:
                model.cploutnam = 'State.Q'
                model.ini_frliq = \
                    model.study.masc.get('State.Z', model.cplsect)

        model.cploutbc = np.zeros(model.tbcinstp, dtype=np.float64)
        model.zrefcpl = model.study.masc.get('Model.Zbot', model.cplsect)
        model.cpl1dcrt = np.zeros(model.nbcriteria, np.double)

        self.mod = model
        self.conv = None

    def init_conv(self):
        """
        Intialisation of convergence variable
        """
        self.conv = ClassConvergence(self.mod.nb_1d_models, self.mod.cplsteps,
                                     output=False)
        self.conv.set_criteria(height=self.mod.crit_arret_h,
                               velocity=self.mod.crit_arret_v)

    def init_cplstp(self):
        """
        Initialise a coupling step
            1. Reset the convergence flag
            2. Update the time bounds of the step
            3. Interpolates the initial BC values on the coupling times
        """

        # Reset the convergence flag
        self.mod.converged = False

        # Update the time bounds of the step
        self.mod.t_i = self.mod.ti_ini + self.mod.cplstp * self.mod.cplfreq
        self.mod.t_f = self.mod.t_i + self.mod.cplfreq
        self.mod.initst = self.mod.cplstp * self.mod.tstinstp
        self.mod.endtst = (self.mod.cplstp + 1) * self.mod.tstinstp

        if self.mod.cpl1dstart == 'persist1d':
            # Recover the last state variables to be used as initial forcings
            if self.mod.cpllaw.type == 1:
                self.mod.cploutnam = 'State.Z'
                self.mod.ini_frliq = self.mod.study.masc.get('State.Q',
                                                             self.mod.cplsect)
            elif self.mod.cpllaw.type == 2:
                self.mod.cploutnam = 'State.Q'
                self.mod.ini_frliq = self.mod.study.masc.get('State.Z',
                                                             self.mod.cplsect)

        # Interpolate the boundary conditions on the current coupling interval
        self.select_bc()

    def select_bc(self):
        """
        Fill the reference boundary conditions array on the coupling interval
            1. Interpolates the initial BC values on the coupling times
        """
        mod = self.mod
        mod.times_bc = [mod.t_i + ib * mod.d_t for ib in range(mod.tbcinstp)]
        mod.times_bc = np.array(mod.times_bc)
        mod.times2d = [mod.t_i + ib * mod.dt_2d for ib in range(mod.nit_2d)]
        mod.times2d = np.array(mod.times2d)
        mod.current_bc = np.zeros((mod.tbcinstp, mod.nbbc), dtype=np.float64)
        for law in range(mod.nbbc):
            mod.current_bc[:, law] = mod.hydlaw[law].time_interp(mod.times_bc)

    def compute_criteria(self):
        """
        Compute the convergence criteria at the coupling sections
         for each 1D model
        Vitesse a l'interface Q/(S1+S2)
        Hauteur d'eau a l'interface Z-ZREF
        Les debits a l'interface Q
        Premier invariant de Riemann Q/(S1+S2)+2*SQRT(G*(Z-ZREF))
        Second invariant de Riemann  Q/(S1+S2)-2*SQRT(G*(Z-ZREF))
        Surface mouillee             S1+S2
        """

        stot = self.mod.study.masc.get(
            'State.S1', self.mod.cplsect) + \
            self.mod.study.masc.get('State.S2', self.mod.cplsect)
        qcpl = self.mod.study.masc.get('State.Q', self.mod.cplsect)
        zcpl = self.mod.study.masc.get('State.Z', self.mod.cplsect)
        self.mod.cpl1dcrt[0] = qcpl / stot
        self.mod.cpl1dcrt[1] = zcpl - self.mod.zrefcpl
        self.mod.cpl1dcrt[2] = qcpl
        self.mod.cpl1dcrt[3] = self.mod.cpl1dcrt[0] + \
            2.0 * sqrt(9.81 * self.mod.cpl1dcrt[1])
        self.mod.cpl1dcrt[4] = self.mod.cpl1dcrt[0] - \
            2.0 * sqrt(9.81 * self.mod.cpl1dcrt[1])
        self.mod.cpl1dcrt[5] = stot

    def update_bc(self):
        """
        Prepare the boundary conditions for the current integration
        1. For the very first integration start from the initial waterline
            values at the coupling sections (constant in time)
        2. For the first iteration on regular coupling steps start
            accordingly to cpl1dstart
           previous2d: start with the last conditions received
            from the 2D model on
           the previous coupling step and interpolate them in time.
           persist2d: start with the last conditions received
            from the 2D model on
           the previous coupling step but keep only the last time step value
            (constant in time)
           persist1d: start with the last 1d state values at the
            coupling sections (constant in time)
        3. For any regular iteration and for the initial iteration at every new
           coupling step, receive the conditions from the 2D model
            and interpolate them in time
        """

        if self.mod.cpl1dstart == 'previous2d':
            self.update_previous2d()
        elif self.mod.cpl1dstart == 'persist2d':
            self.update_persist2d()
        elif self.mod.cpl1dstart == 'persist1d':
            self.update_persist1d()

    #
    #        Optionally print out the updated boundary condition
    #        if self.mod.config == 'NA':
    #            print ('Stp ',self.mod.cplstp,' It ',self.mod.iter,' Bc ', \
    #                       self.mod.current_bc[:,self.mod.cpllaw.pos])
    #

    def update_persist1d(self):
        """
        Update boundary conditions with persist1d method:
        start with the last 1d state values at the coupling
        sections (constant in time)
        """
        if self.mod.iter == 0:
            if self.mod.cplstp > 0:
                self.mod.comm.receive_bc()
            self.mod.current_bc[:, self.mod.cpllaw.pos] = \
                copy.deepcopy(self.mod.ini_frliq)
        else:
            self.mod.comm.receive_bc()
            self.mod.current_bc[:, self.mod.cpllaw.pos] = np.interp(
                self.mod.times_bc, self.mod.times2d, self.mod.cpl_from2d)

    def update_persist2d(self):
        """
        Update boundary conditions with persist2d method:
        Start with the last conditions received from the 2D model
         on the previous
        coupling step but keep only the last time step value (constant in time)

        """
        if self.mod.cplstp == 0 and self.mod.iter == 0:
            if self.mod.inirun:
                self.mod.current_bc[:, self.mod.cpllaw.pos] = \
                    self.mod.ini_frliq
                self.mod.last_from2d = copy.deepcopy(self.mod.ini_frliq)
            else:
                if os.path.exists(self.mod.jsonbc):
                    with open(self.mod.jsonbc) as json_data:
                        dico = json.load(json_data)
                        try:
                            self.mod.last_from2d = \
                                copy.deepcopy(dico['last_from2d'])
                        except KeyError:
                            self.mod.last_from2d = \
                                copy.deepcopy(dico['conlim2d'][-1])
                        self.mod.last_from2d = float(self.mod.last_from2d)
                        self.mod.cpl_from2d = \
                            np.array(dico['conlim2d'], dtype=np.float64)
                    self.mod.last_from2d = \
                        (self.mod.cpl_from2d[-1] + self.mod.last_from2d) * 0.5
                    self.mod.current_bc[:, self.mod.cpllaw.pos] = \
                        self.mod.last_from2d
                    os.remove(self.mod.jsonbc)
                else:
                    self.mod.current_bc[:, self.mod.cpllaw.pos] = \
                        self.mod.ini_frliq
                    self.mod.last_from2d = copy.deepcopy(self.mod.ini_frliq)
        elif self.mod.cplstp > 0 and self.mod.iter == 0:
            self.mod.comm.receive_bc()
            # # Restart from the last time instance of the 2d law
            # self.mod.current_bc[:,self.mod.cpllaw.pos] =
            #   self.mod.cpl_from2d[-1]
            # # Restart from the time average of the
            # 2d law on the previous cpl step
            # self.mod.current_bc[:,self.mod.cpllaw.pos] =
            #   np.average(self.mod.cpl_from2d)
            # Restart from the average of the values of
            # the last time instance of
            # the 2d law from the two last iterates on the previous cpl stp
            self.mod.last_from2d = (self.mod.cpl_from2d[-1] +
                                    self.mod.last_from2d) * 0.5
            self.mod.current_bc[:, self.mod.cpllaw.pos] = self.mod.last_from2d
        else:
            self.mod.comm.receive_bc()
            self.mod.last_from2d = copy.deepcopy(self.mod.cpl_from2d[-1])
            self.mod.current_bc[:, self.mod.cpllaw.pos] = np.interp(
                self.mod.times_bc, self.mod.times2d, self.mod.cpl_from2d)

    def update_previous2d(self):
        """
        Update boundary conditions with previous2d method:
        start with the last conditions received from the 2D model on
        the previous coupling step and interpolate them in time.
        """
        if self.mod.cplstp == 0 and self.mod.iter == 0:
            if self.mod.inirun:
                self.mod.current_bc[:, self.mod.cpllaw.pos] = \
                    copy.deepcopy(self.mod.ini_frliq)
            else:
                if os.path.exists(self.mod.jsonbc):
                    with open(self.mod.jsonbc) as json_data:
                        dico = json.load(json_data)
                        self.mod.cpl_from2d = np.array(dico['conlim2d'],
                                                       dtype=np.float64)
                    self.mod.current_bc[:, self.mod.cpllaw.pos] = np.interp(
                        self.mod.times_bc, self.mod.times2d,
                        self.mod.cpl_from2d)
                    os.remove(self.mod.jsonbc)
                else:
                    self.mod.current_bc[:, self.mod.cpllaw.pos] = \
                        copy.deepcopy(self.mod.ini_frliq)
        else:
            self.mod.comm.receive_bc()
            self.mod.current_bc[:, self.mod.cpllaw.pos] = np.interp(
                self.mod.times_bc, self.mod.times2d, self.mod.cpl_from2d)

    def transmit_bc(self):
        """
        Transmit the computed boundary conditions to the 2D model
        """
        # Transmit BC_1D to 2d
        self.mod.comm.transmit_bc()

    def store_bc(self):
        """
        Store the 1D BC (current and previous) for restarting
        """
        # Transmit BC_1D to 2d
        self.mod.comm.receive_bc()

        if self.mod.cpl1dstart == 'persist2d':
            dico = {'model': self.mod.config,
                    'last_from2d': self.mod.last_from2d,
                    'conlim2d': self.mod.cpl_from2d.tolist()}
        else:
            dico = {'model': self.mod.config,
                    'conlim2d': self.mod.cpl_from2d.tolist()}
        jsonbcout = 'bc1D_restart_' + self.mod.config + '_' + \
            str(int(self.mod.t_f)) + '.json'
        with open(jsonbcout, 'w') as outfile:
            json.dump(dico, outfile)

    def check_convergence(self):
        """
        Check if the current iteration has reached convergence
        """

        # Transmit the 1D criteria
        self.mod.comm.transmit_criteria()
        # Get the 2D criteria
        self.mod.comm.receive_criteria()
        # Compute convergence flag on root
        if self.mod.ismaster:
            self.mod.converged = \
                self.conv.main(self.mod.iter,
                               self.mod.vars_1d, self.mod.vars_2d)
        # Synchro of convergence flags
        self.mod.comm.bcast_convergence()


class HydLaw:
    """
    hydraulic law processing
    """

    def __init__(self, model, pos, name, idxhyd):
        """
        Constructor
            1. Initialize empty values
            2. Recover the initial time and values as from
                the prescribed chronicles
        @param model (class_mod_1d) pointer to 1D model instance
        @param pos (int) position in the BC list
        @param name (str) name of the associated hydraulic law
        @param idxhyd (int) number of the associated hydraulic law
        """

        self.pos = pos
        self.name = name
        self.num = idxhyd
        self.size = 0

        self.type = model.study.masc.get('Model.Graph.Type', self.num)
        if self.type == 1:
            varname = 'Model.Graph.Discharge'
        elif self.type == 2:
            varname = 'Model.Graph.Level'

        _, self.size, _ = model.study.masc.get_var_size(varname, self.num)
        self.bc_ini = np.ones(self.size, dtype=np.float64)
        self.bc_time = np.ones(self.size, dtype=np.float64)
        for idt in range(self.size):
            self.bc_ini[idt] = model.study.masc.get(varname, self.num, idt)
            self.bc_time[idt] = \
                model.study.masc.get('Model.Graph.Time', self.num, idt)

    def time_interp(self, times):
        """
        Provides the values of the law at the prescribed coupling times
            1. Uses the numpy interpolator
        @param times (ndarray) numpy array of prescribed times
        @return (np.array) new law
        """

        #  newlaw = np.zeros(times.size, dtype=np.float64)
        newlaw = np.interp(times, self.bc_time, self.bc_ini)

        return newlaw
