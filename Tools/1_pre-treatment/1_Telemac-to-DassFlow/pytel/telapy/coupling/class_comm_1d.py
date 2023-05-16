# coding: utf-8
"""
ClassComm1D class
=================
"""
import numpy as np
from mpi4py import MPI


class ClassComm1D:
    """
    Communication tools for the 1D models
        1. Split the global MPI communicator into its 1D and 2D subcomm
        2. Compute and store the rank in the communicator
    """

    def __init__(self, model):
        """
        Constructor
        @param model (obj) : ClassMod1D model, pointer to 1D model instance
        """
        self.cpl_comm = MPI.COMM_WORLD
        self.comm_1d = self.cpl_comm.Split(1)
        self.rank = self.comm_1d.Get_rank()

        # Associate model instance
        self.mod = model
        self.mod.instance = self.rank  # True if one 1D model per process
        self.mod.ismaster = self.rank == 0

        # communication variables
        self.bc_1to2 = None
        self.bc_2to1 = None
        self.cr_1to2 = None
        self.cr_2to1 = None

    def init_comms(self):
        """
        Initialize persistent MPI communications for recurrent communications.
        """

        if self.rank == 0:
            # Transmitted boundary conditions
            self.mod.bcsendvar = \
                np.zeros([self.mod.nb_1d_models, self.mod.tbcinstp],
                         dtype=np.float64)
            dest = self.comm_1d.Get_size()
            # 2D has the first rank after 1D instances
            self.bc_1to2 = \
                self.cpl_comm.Send_init(self.mod.bcsendvar, dest, 1012)

            # Received boundary conditions
            self.mod.bcgetvar = \
                np.zeros([self.mod.nb_1d_models, self.mod.nit_2d],
                         dtype=np.float64)
            ori = self.comm_1d.Get_size()
            # 2D has the first rank after 1D instances
            self.bc_2to1 = \
                self.cpl_comm.Recv_init(self.mod.bcgetvar, ori, 1021)

            # Transmitted 1D convergence criteria
            self.mod.vars_1d = \
                np.zeros([self.mod.nb_1d_models, self.mod.nbcriteria],
                         dtype=np.float64)
            dest = self.comm_1d.Get_size()
            # 2D has the first rank after 1D instances
            self.cr_1to2 = \
                self.cpl_comm.Send_init(self.mod.vars_1d, dest, 2012)

            # Received 2D convergence criteria
            self.mod.vars_2d = \
                np.zeros([self.mod.nb_1d_models, self.mod.nbcriteria],
                         dtype=np.float64)
            ori = self.comm_1d.Get_size()
            # 2D has the first rank after 1D instances
            self.cr_2to1 = self.cpl_comm.Recv_init(self.mod.vars_2d, ori, 2021)

    def receive_bc(self):
        """
        Receives the boundary conditions at the interfaces from the 2D model.

        1. Receive the whole set
        2. Scatter the components to the individual models
        """

        self.mod.cpl_from2d = np.zeros(self.mod.nit_2d, dtype=np.float64)
        if self.rank == 0:
            self.bc_2to1.Start()
            self.bc_2to1.wait()
            self.comm_1d.Scatter(self.mod.bcgetvar, self.mod.cpl_from2d,
                                 root=0)
        else:
            self.comm_1d.Scatter(None, self.mod.cpl_from2d, root=0)

    def transmit_bc(self):
        """
        Transmit the computed boundary conditions to the 2D model
        1. Collect the computed BC's on the master processor
        2. Transmit them
        """

        # Collect on master proc the 1d bc's of all 1d models and
        # transmit BC_1D to 2d
        if self.rank == 0:
            self.comm_1d.Gather(self.mod.cploutbc, self.mod.bcsendvar, root=0)
            self.bc_1to2.Start()
            self.bc_1to2.Wait()
        else:
            self.comm_1d.Gather(self.mod.cploutbc, None, root=0)

    def transmit_criteria(self):
        """
        Transmit the computed criteria at coupling section
         to the convergence unit

        1. Collect the computed criteria on the master processor
        2. Transmit them
        """

        # Collect on master proc the 1d bc's of all 1d models and
        # transmit criteria
        if self.rank == 0:
            self.comm_1d.Gather(self.mod.cpl1dcrt, self.mod.vars_1d, root=0)
            self.cr_1to2.Start()
            self.cr_1to2.Wait()
        else:
            self.comm_1d.Gather(self.mod.cpl1dcrt, None, root=0)

    def receive_criteria(self):
        """
        Receive the criteria variables
        """
        if self.rank == 0:
            self.cr_2to1.Start()
            self.cr_2to1.Wait()

    def bcast_convergence(self):
        """
        Function allows to Bcast the convergence variables
        """
        self.mod.converged = self.comm_1d.bcast(self.mod.converged, root=0)
