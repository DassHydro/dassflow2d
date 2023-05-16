# -*- coding: utf-8 -*-
"""
ClassCpl2D
"""
from math import sqrt

import numpy as np
from mpi4py import MPI


class ClassCpl2D:
    """
    Coupling tools for the 2D models
    """

    def __init__(self, t2d, coupler, nb_model_1d, nb_criteria, id_fr_co_2d):
        """
        Constructor:
        Initialization of Telemac instance
        Initialization of coupler variables
        Initialization of coupling model variables

        @param t2d (Telemac2d) Telemac2d instance
        @param coupler (???) Coupler structure
        @param nb_model_1d (int) Number of 1d models
        @param nb_criteria (int) Number of criteria
        @param id_fr_co_2d (np.array) border point index
        """
        self.nb_model_1d = nb_model_1d
        self.nb_criteria = nb_criteria
        self.coupler = coupler
        self.comm = self.coupler.comm_2d
        self.t2d = t2d
        # initialise the coupling interface descriptors
        self.node_frliq(id_fr_co_2d)

        self.glo_zf = self.gather_fr_co('MODEL.BOTTOMELEVATION')
        self.glo_x = self.gather_fr_co('MODEL.X')
        self.glo_y = self.gather_fr_co('MODEL.Y')

        self.glo_h = list()
        self.flux_boundaries = list()
        self.interp_conlim_co = list()

        if self.coupler.rank == 0:
            self.mini_zf = []
            for i in range(self.nb_fr_co):
                self.mini_zf.append(np.min(self.glo_zf[i]))
            self.mini_zf = np.array(self.mini_zf)

    def node_frliq(self, id_fr_co_2d):
        """
        Search frontiere  liquide
        @param id_fr_co_2d (np.array) border point index
        """
        self.id_fr_co_2d = id_fr_co_2d.astype(int)
        self.nb_fr_co = len(self.id_fr_co_2d)
        self.tab_frliq = []
        self.nb_node_loc = np.zeros(self.nb_fr_co, dtype=np.int)
        self.nb_node_glo = np.zeros(self.nb_fr_co, dtype=np.int)
        nptfr = self.t2d.get('MODEL.NPTFR')
        nbor_numliq = self.t2d.get_array('MODEL.NBOR')
        numliq = self.t2d.get_array('MODEL.NUMLIQ')
        numliq = numliq.astype(int)
        # Always include the comm_2d master in the restricted
        # coupling interfaces comm
        self.has_fr_co = self.coupler.rank == 0
        # Mesh internal boundaries processor association
        nachb = {}
        if self.coupler.ncsize > 1:
            nptir = self.t2d.get('MODEL.NPTIR')
            nachb_l = self.t2d.get_array('MODEL.NACHB')
            for i in range(nptir):
                nachb[nachb_l[i, 0]] = nachb_l[i, 1] - 1
        for i in range(self.nb_fr_co):
            self.get_node_frco(i, nptfr, numliq, nbor_numliq)
            self.del_duplic_halo_pt(i, nachb)

            self.has_fr_co = self.has_fr_co or len(self.tab_frliq[i]) > 0
            self.nb_node_loc[i] = len(self.tab_frliq[i])
        if self.has_fr_co:
            self.tab_frliq = np.array(self.tab_frliq)
        # create a small communicator for the processors
        # having some coupling points
        self.coupler.comm_fr_co = self.comm.Split(self.has_fr_co)
        self.comm_fr_co = self.coupler.comm_fr_co

        # the master needs to kwnow how many points each coupling interface has
        for i in range(self.nb_fr_co):
            if self.coupler.ncsize > 1:
                if self.has_fr_co:
                    val_tmp = self.comm_fr_co.reduce(self.nb_node_loc[i],
                                                     op=MPI.SUM, root=0)
                    if self.coupler.rank == 0:
                        self.nb_node_glo[i] = val_tmp
            else:
                self.nb_node_glo[i] = len(self.tab_frliq[i])

    def del_duplic_halo_pt(self, idx, nachb):
        """
        Detect duplicated halo points and delete that
        @param idx (int) iteration step
        @param nachb (dict) Mesh internal boundaries processor association
        """
        # detect duplicated halo points
        if self.coupler.ncsize > 1:
            for k in reversed(range(len(self.tab_frliq[idx]))):
                if self.tab_frliq[idx][k] in nachb:
                    if self.coupler.rank != nachb[self.tab_frliq[idx][k]]:
                        del self.tab_frliq[idx][k]

    def get_node_frco(self, idx, nptfr, numliq, nbor_numliq):
        """
        Get nodes of coupling liquid border (self.tab_frliq)
        @param idx (int) iteration step
        @param nptfr (int) liquid border number
        @param numliq (list) list of liquid border number
        @param nbor_numliq (list) list of liquid border node
        """
        compteur_g = 0
        tab_node_frliq_tmp = []
        tmp2 = []
        # nb of border point on that straddling BC
        if nptfr > 0:
            if (numliq[0] == self.id_fr_co_2d[idx]) and \
                    (numliq[nptfr - 1] == self.id_fr_co_2d[idx]):
                while numliq[nptfr - compteur_g - 1] != 0:
                    compteur_g += 1

        for j in range(nptfr):
            # if J point is before the straddling border
            if (j + compteur_g) <= nptfr - 1:
                if numliq[j] == self.id_fr_co_2d[idx]:
                    tab_node_frliq_tmp.append(int(nbor_numliq[j]))
            # if J point is along the straddling border
            elif (j + compteur_g) > nptfr - 1:
                if numliq[j] == self.id_fr_co_2d[idx]:
                    tmp2.append(int(nbor_numliq[j]))
        self.tab_frliq.append(tmp2 + tab_node_frliq_tmp)

    @staticmethod
    def compt_type_cl(type_cl_in_1d):
        """
        Compute the number of boundary conditions type
        @param type_cl_in_1d (list) type of boundary condition
        @return ncote (int) nb of  level water condition
        @return ndeb (int) nb of flow rate condition
        """
        unique, counts = np.unique(type_cl_in_1d, return_counts=True)
        tmp = {x: y for x, y in zip(unique, counts)}
        if 1 not in tmp.keys():
            tmp[1] = 0
        if 2 not in tmp.keys():
            tmp[2] = 0

        return tmp[1], tmp[2]

    def interp_conlim(self, nit, d_t, dt_1d, conlim_co):
        """
        Temporal interpolation of boundary condition
        @param nit (int) iteration number of 2D model
        @param d_t (float) Time step of 2D model
        @param dt_1d (float) Time step of 1D model
        @param conlim_co (np.array) boundary condition table
        """
        if self.has_fr_co:
            self.interp_conlim_co = np.zeros([nit + 1, self.nb_model_1d])
            tps2d = [v * d_t for v in range(nit + 1)]
            tps1d = [v * dt_1d for v in range(len(conlim_co[0, :]))]

            for j in range(self.nb_model_1d):
                self.interp_conlim_co[:, j] = \
                    np.interp(tps2d, tps1d, conlim_co[j, :])

    def comput_conlim(self, pos_model_1d, type_cl_in_1d):
        """
        Compute boudaries conditions
        @param type_cl_in_1d (list) type of boundary conditions
        @param pos_model_1d (list) models 1D position
        @return conlim (np.array) boundary condition values
        """
        if 2 in type_cl_in_1d:
            self.glo_h = self.gather_fr_co('MODEL.WATERDEPTH')

        if self.coupler.rank == 0:
            if 1 in type_cl_in_1d:
                self.flux_boundaries = \
                    self.t2d.get_array('MODEL.FLUX_BOUNDARIES')
            conlim = []
            for i, typec in enumerate(type_cl_in_1d):
                if typec == 1:
                    if pos_model_1d[i] == 2:
                        # id_fr_co_2d[i]-1 because of  fortran begin 1
                        conlim.append(self.flux_boundaries[self.id_fr_co_2d[i]
                                      - 1])
                    elif pos_model_1d[i] == 1:
                        conlim.append(-self.flux_boundaries[self.id_fr_co_2d[i]
                                      - 1])
                elif typec == 2:
                    conlim.append(self.calc_sl_moy(self.glo_h, i))
        else:
            conlim = None

        return conlim

    def calc_sl_moy(self, glo_h, i=0):
        """
        Compute mean free surface
        @param glo_h (np.array) height values of global mesh
        @param i (int) numero of 1D model
        @return (float) Mean free surface
        """
        moy = []
        for k in range(1, self.nb_node_glo[i] - 1):
            if glo_h[i][k] > 1e-3:
                moy.append(glo_h[i][k] + self.glo_zf[i][k])
        return np.mean(np.array(moy))

    def cl_modif(self, type_cl_in_1d, condlim):
        """
        update boundaries conditions
        @param type_cl_in_1d (list) type of boundary condition
        @param condlim (np.array) boundary conditions
        """
        if self.has_fr_co:
            cote = self.t2d.get_array('MODEL.COTE')
            deb = self.t2d.get_array('MODEL.DEBIT')
            for i, typc in enumerate(type_cl_in_1d):
                # ! cas 1 : The 1D model have Q in Boundary
                # condition(BC), it send  Z
                if typc == 1:
                    cote[self.id_fr_co_2d[i] - 1] = condlim[i]
                    # ! cas 2 :  The 1D model have Z in BC, it send  Q
                elif typc == 2:
                    deb[self.id_fr_co_2d[i] - 1] = condlim[i]
            self.t2d.set_array('MODEL.COTE', cote)
            self.t2d.set_array('MODEL.DEBIT', deb)

    def calc_area(self, x_a, x_b, z_m, epsi=0.001):
        """
        Compute area of one mesh

        @param x_a (list) [ x,y position for A point, bottom for A point]
        @param x_b (list)  [ x,y position for B point, bottom for B point]
        @param z_m (float)  Mean height
        @param epsi (float) tolerance by default it's equal 0.001

        @return air (float) return air
        """
        air = 0
        if z_m - x_a[1] >= epsi and z_m - x_b[1] >= epsi:
            air = (2 * z_m - x_a[1] - x_b[1]) * self.dist(x_b[0], x_a[0]) * 0.5
            # It is considered that the water level is straight on the elements
            # at the ends of the liquid border. What is not
            # the case in Telemac.
        elif z_m - x_a[1] < epsi <= z_m - x_b[1]:
            air = self.dist(x_b[0], x_a[0]) * (z_m - x_b[1]) ** 2 / \
                (2 * (x_a[1] - x_b[1]))
        elif z_m - x_a[1] >= epsi > z_m - x_b[1]:
            air = self.dist(x_b[0], x_a[0]) * (z_m - x_a[1]) ** 2 / \
                (2 * (x_b[1] - x_a[1]))
        return air

    @staticmethod
    def dist(x_a, x_b):
        """
        Compute distance

        @param x_a (list) a point coord.
        @param x_b (list) b point coord.

        @return (float) the distance between a point and b point
        """
        return sqrt((x_a[0] - x_b[0]) ** 2 + (x_a[1] - x_b[1]) ** 2)

    def var_interface(self, pos_model_1d, type_cl_in_1d, vars_2d):
        """
        Creating the table to calculate the convergence variable

        @param pos_model_1d (list) models 1D position
        @param type_cl_in_1d (list) boundary condition type
        @param vars_2d (np.array) table to calculate the convergence variable
        """

        if 2 not in type_cl_in_1d:
            self.glo_h = self.gather_fr_co('MODEL.WATERDEPTH')

        if self.coupler.rank == 0:

            if 1 not in type_cl_in_1d:
                self.flux_boundaries = \
                    self.t2d.get_array('MODEL.FLUX_BOUNDARIES')

            grav = 9.81

            for i, pos in enumerate(pos_model_1d):
                # mean surface
                z_2d = self.calc_sl_moy(self.glo_h, i)
                area_sm = 0
                for idb in range(self.nb_node_glo[i] - 1):
                    idb2 = idb + 1
                    # we consider THE FREE SURFACE ON THE LIQUID BORDER IS
                    # straight to remain consistent with Mascaret
                    x_1 = [self.glo_x[i][idb], self.glo_y[i][idb]]
                    x_2 = [self.glo_x[i][idb2], self.glo_y[i][idb2]]
                    area_sm += self.calc_area(
                        [x_1, self.glo_zf[i][idb]],
                        [x_2, self.glo_zf[i][idb2]],
                        z_2d)
                # calcul flux
                v_2d = None
                if pos == 2:
                    # id_fr_co_2d[i]-1 because of  fortran begin 1
                    vars_2d[i, 2] = self.flux_boundaries[self.id_fr_co_2d[i]
                                                         - 1]
                    v_2d = vars_2d[i, 2] / area_sm
                elif pos == 1:
                    vars_2d[i, 2] = -self.flux_boundaries[self.id_fr_co_2d[i]
                                                          - 1]
                    v_2d = vars_2d[i, 2] / area_sm

                # Velocity into interface
                vars_2d[i, 0] = v_2d

                # Water height into interface
                vars_2d[i, 1] = z_2d - self.mini_zf[i]
                # compute Riemann invrarient
                inva_i2d = vars_2d[i, 0] + 2 * sqrt(grav * vars_2d[i, 1])
                inva_j2d = vars_2d[i, 0] - 2 * sqrt(grav * vars_2d[i, 1])
                vars_2d[i, 3] = inva_i2d
                vars_2d[i, 4] = inva_j2d
                vars_2d[i, 5] = area_sm
        else:
            vars_2d = None
        return vars_2d

    def gather_fr_co(self, vars_n):
        """
        Gather boundary conditions

        @param vars_n (str)  Variable name
        """

        if self.has_fr_co:
            loc2d = self.t2d.get_array(vars_n)
            if self.coupler.rank == 0:
                valgr = []
            else:
                valgr = None
            for i in range(self.nb_fr_co):
                valtmp = [loc2d[self.tab_frliq[i][j] - 1]
                          for j in range(self.nb_node_loc[i])]
                if self.coupler.ncsize > 1:
                    vallo = self.comm_fr_co.reduce(valtmp, op=MPI.SUM, root=0)
                else:
                    vallo = valtmp
                if self.coupler.rank == 0:
                    valgr.append(vallo)
        else:
            valgr = None

        return valgr
