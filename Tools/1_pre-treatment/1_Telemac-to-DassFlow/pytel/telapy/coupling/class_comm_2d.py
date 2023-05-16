# -*- coding: utf-8 -*-
"""
ClassComm2D class
=================
"""
import numpy as np
from mpi4py import MPI


class ClassComm2D:
    """
    Communication tools for the 2D models
        1. Split the global MPI communicator into its 1D and 2D subcomm
        2. Compute and store the rank in the communicator
    """

    def __init__(self, model):
        """
        Constructor
        @param model (obj) : ClassMod2D model, pointer to 2D model instance
        """
        self.cpl_comm = MPI.COMM_WORLD
        self.comm_2d = self.cpl_comm.Split(2)
        # Useful only when T2D will be parallel
        self.rank = self.comm_2d.Get_rank()
        self.ncsize = self.comm_2d.Get_size()
        self.mod = model

        self.bc_1to2 = None
        self.bc_2to1 = None
        self.cr_1to2 = None
        self.cr_2to1 = None

        self.comm_fr_co = None

    def init_comms(self):
        """
        Initialize persistent MPI communications for recurrent communications.
        """
        self.mod.conlim_2d = np.zeros([self.mod.nb_model_1d, self.mod.nit + 1],
                                      dtype=np.float64)
        self.mod.conlim_co = np.zeros([self.mod.nb_model_1d, self.mod.nit_1d],
                                      dtype=np.float64)
        # Transmitted 2D convergence criteria
        self.mod.vars_2d = \
            np.zeros([self.mod.nb_model_1d, self.mod.nb_criteria],
                     dtype=np.float64)
        # Received 1D convergence criteria
        self.mod.vars_1d = \
            np.zeros([self.mod.nb_model_1d, self.mod.nb_criteria],
                     dtype=np.float64)

        if self.rank == 0:  # Useful only when T2D will be parallel
            # Transmitted boundary conditions
            dest = 0  # 1D instances are before 2D in mpirun
            self.bc_2to1 = \
                self.cpl_comm.Send_init(self.mod.conlim_2d, dest, 1021)

            # Received boundary conditions
            ori = 0  # 1D instances are before 2D in mpirun
            self.bc_1to2 = \
                self.cpl_comm.Recv_init(self.mod.conlim_co, ori, 1012)

            # Transmitted 2D convergence criteria
            dest = 0  # 1D instances are before 2D in mpirun
            self.cr_2to1 = \
                self.cpl_comm.Send_init(self.mod.vars_2d, dest, 2021)

            # Received 1D convergence criteria
            ori = 0  # 1D instances are before 2D in mpirun
            self.cr_1to2 = self.cpl_comm.Recv_init(self.mod.vars_1d, ori, 2012)

    def get_cl(self):
        """
        get boundaries conditions values
        """
        if self.rank == 0:  # Useful only when T2D will be parallel
            self.bc_1to2.Start()
            self.bc_1to2.Wait()
        if self.ncsize > 1:
            conlim_tmp = self.comm_fr_co.bcast(self.mod.conlim_co, root=0)
            if self.rank > 0:
                self.mod.conlim_co = np.array(conlim_tmp)

    def transmit_cl_1d(self):
        """
        Push conlim_2d and vars_2d variables
        """
        # send
        if self.rank == 0:  # Useful only when T2D will be parallel
            # send criteria
            self.cr_2to1.Start()
            self.cr_2to1.Wait()
            # send BC variables
            self.bc_2to1.Start()
            self.bc_2to1.Wait()

    def get_vars_1d(self):
        """
        Get 1D criteria  variables
        """
        if self.rank == 0:  # Useful only when T2D will be parallel
            self.cr_1to2.Start()
            self.cr_1to2.Wait()
