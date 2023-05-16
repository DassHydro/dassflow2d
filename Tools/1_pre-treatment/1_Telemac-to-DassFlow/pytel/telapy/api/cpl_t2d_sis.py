#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
    Class for a t2d sis coupled model

    Author(s): Fabrice Zaoui, Yoann Audouin, Cedric Goeury, Renaud Barate

    Copyright EDF 2016
"""
from telapy.api.t2d import Telemac2d
from telapy.api.sis import Sisyphe


class CplT2dSis(object):
    """The Generic Python class for t2d-sis coupling"""

    def __init__(self, t2d_steering_file,
                 sis_steering_file,
                 user_fortran, stdout=6,
                 log_lvl='INFO',
                 comm=None, recompile=False):
        """
        Constructor for apiModule

        @param t2d_steering_file (string) Name of the Telemac2d steering file
        @param sis_steering_file (string) Name of the Sisyphe steering file
        @param user_fortran (string) Name of the user Fortran
        @param stdout (int) Where to put the listing
        @param comm (MPI.Comm) MPI communicator
        @param recompile (boolean) If true recompiling the API
        @param log_lvl (string) Logger level
        """
        self.t2d = Telemac2d(t2d_steering_file,
                             user_fortran=user_fortran,
                             recompile=True,
                             stdout=stdout,
                             log_lvl=log_lvl,
                             comm=comm)

        self.sis = Sisyphe(sis_steering_file,
                           stdout=stdout,
                           comm=comm,
                           log_lvl=log_lvl,
                           recompile=False)

    def set_case(self):
        """
           Read the steering file and run allocation
        """
        self.t2d.set_case()
        self.sis.set_case(init=False)

    def init_state_default(self):
        """
        Initialize the state of the model Telemac 2D with the values of
        disharges and water levels as indicated by the steering file
        """
        # Init of Telemac2d
        self.t2d.init_state_default()
        # Init of coupling
        self.t2d.api_inter.cpl_init(self.t2d.my_id)
        # transfering sis data
        self.t2d.api_inter.set_var_sis(self.t2d.my_id, self.sis.my_id, 0)
        # Init of sisyphe
        self.sis.init_state_default()
        # Copying back data to telemac2d
        self.t2d.api_inter.set_var_t2d(self.t2d.my_id, self.sis.my_id)
        # Saving charr_susp state
        self.t2d.api_inter.save_charr_susp(self.t2d.my_id, self.sis.my_id)
        #
        self.cpl_period = self.t2d.get('MODEL.CPL_PERIOD')

    def run_one_time_step(self, istep):
        """
        Run one time step
        """
        self.t2d.api_inter.run_timestep_compute_t2d(self.t2d.my_id)
        if(self.cpl_period*(int(istep/self.cpl_period)) == istep):
            self.t2d.api_inter.run_timestep_sis_cpl(self.t2d.my_id,
                                                    self.sis.my_id)
        self.t2d.api_inter.run_timestep_res_t2d(self.t2d.my_id)

    def run_all_time_steps(self):
        """
        Run all the time steps

        @return the number of computed time steps
        """

        ntimesteps = self.t2d.get('MODEL.NTIMESTEPS')
        for istep in range(ntimesteps):
            self.run_one_time_step(istep)

    def finalize(self):
        """
        Delete the Telemac 2D instance

        @return error code
        """
        self.sis.api_inter.run_finalize_sis(self.sis.my_id)
        self.t2d.api_inter.run_finalize_t2d(self.t2d.my_id)

    def __del__(self):
        del self.sis
        del self.t2d
