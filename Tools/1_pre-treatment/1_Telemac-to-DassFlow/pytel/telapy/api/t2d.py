# -*- coding: utf-8 -*-
"""
    Python wrapper to the Fortran APIs of Telemac 2D

    Author(s): Fabrice Zaoui, Yoann Audouin, Cedric Goeury, Renaud Barate

    Copyright EDF 2016
"""

import os
import numpy as np
from telapy.api.api_module import ApiModule
from utils.exceptions import TelemacException
from utils.polygon import is_in_polygon


class Telemac2d(ApiModule):
    """The Telemac 2D Python class for APIs"""
    _instanciated = False

    def __new__(cls, *args, **kwargs):
        if cls._instanciated:
            raise TelemacException("a Telemac2d instance already exists")
        instance = ApiModule.__new__(cls)
        cls._instanciated = True
        return instance

    def __init__(self, casfile,
                 user_fortran=None,
                 dicofile=None,
                 lang=2, stdout=6,
                 comm=None,
                 log_lvl='INFO',
                 recompile=True):
        """
        Constructor for Telemac2d

        @param casfile (string) Name of the steering file
        @param user_fortran (string) Name of the user Fortran
        @param dicofile (string) Path to the dictionary
        @param lang (int) Language for ouput (1: French, 2:English)
        @param stdout (int) Where to put the listing
        @param comm (MPI.Comm) MPI communicator
        @param recompile (boolean) If true recompiling the API
        @param log_lvl (string) Logger level
        """
        if dicofile is None:
            hometel = os.getenv("HOMETEL")
            if hometel is not None:
                default_dicofile = os.path.join(os.getenv("HOMETEL"),
                                                "sources",
                                                "telemac2d",
                                                "telemac2d.dico")
            else:
                default_dicofile = 'telemac2d.dico'

            dicofile = default_dicofile
        super(Telemac2d, self).__init__("t2d", casfile, user_fortran,
                                        dicofile, lang, stdout,
                                        comm, recompile, log_lvl=log_lvl)
        self.hsave = None
        self.usave = None
        self.vsave = None
        self.depth = None
        self.u_vel = None
        self.v_vel = None

        self.itrub = None
        self.nbTrac = None
        self.aksave = None
        self.epsave = None
        self.tracsave = None

    def save_state(self):
        """
        Save the hydraulic state
        """
        self.itrub = self.get('MODEL.ITURB')
        self.nbTrac = self.get('MODEL.NTRAC')
        self.hsave = self.get_array('MODEL.WATERDEPTH')
        self.usave = self.get_array('MODEL.VELOCITYU')
        self.vsave = self.get_array('MODEL.VELOCITYV')
        if self.itrub == 3:
            self.aksave = self.get_array('MODEL.AK')
            self.epsave = self.get_array('MODEL.EP')

        if self.nbTrac > 0:
            for k in range(self.nbTrac):
                self.tracsave[k] = self.get_array('MODEL.TRACER')

    def restore_state(self):
        """
        Restore the hydraulic state
        """
        if self.hsave is None:
            raise Exception('Error: unable to restore the hydraulic state.'
                            '\nNo saved state found')
        self.set_array('MODEL.WATERDEPTH', self.hsave)
        self.set_array('MODEL.VELOCITYU', self.usave)
        self.set_array('MODEL.VELOCITYV', self.vsave)
        if self.itrub == 3:
            self.set_array('MODEL.AK', self.aksave)
            self.set_array('MODEL.EP', self.epsave)
        if self.nbTrac > 0:
            for k in range(self.nbTrac):
                self.set_array('MODEL.TRACER', self.tracsave[k])

    def get_state(self):
        """
        Get the hydraulic state

        @returns the hydraulic state: depth (m) .. u_vel (m/s) .. v_vel (m/s)
        """
        self.depth = self.get_array('MODEL.WATERDEPTH')
        self.u_vel = self.get_array('MODEL.VELOCITYU')
        self.v_vel = self.get_array('MODEL.VELOCITYV')
        return self.depth, self.u_vel, self.v_vel

    def set_state(self, hval, uval, vval):
        """
        Set the hydraulic state: hval (m) .. uval (m/s) .. vval (m/s)

        @param hval Water depth value
        @param uval Velocity U value
        @param vval Velocity V value
        """
        self.set_array('MODEL.WATERDEPTH', hval)
        self.set_array('MODEL.VELOCITYU', uval)
        self.set_array('MODEL.VELOCITYV', vval)

    def show_state(self, show=True):
        """
        Show the hydraulic state with matplotlib

        @param show Display the graph (Default True)

        @returns the figure object
        """
        import matplotlib.pyplot as plt
        if self.coordx is not None:
            _, _, _ = self.get_mesh()
        values = self.get_state()
        fig = plt.figure()
        plt.subplot(1, 2, 1)  # water levels
        plt.tripcolor(self.coordx, self.coordy, self.tri, values[0],
                      shading='gouraud', cmap=plt.cm.winter)
        plt.colorbar()
        plt.title('Water levels (m)')
        plt.xlabel('X-coordinate (m)')
        plt.ylabel('Y-coordinate (m)')
        plt.subplot(1, 2, 2)  # velocity
        uvnorm = np.sqrt(values[1]**2 + values[2]**2)
        plt.quiver(self.coordx, self.coordy, values[1], values[2], uvnorm,
                   units='xy', angles='uv', scale=0.01)
        plt.colorbar()
        plt.title('Velocity field (m/s)')
        plt.xlabel('X-coordinate (m)')
        plt.ylabel('Y-coordinate (m)')
        if show:
            plt.show()
        return fig

    def set_bathy(self, bathy, polygon=None):
        """
        Set a new bathy in the geometry file

        @param bathy Array containing the new bathymetry for each point
        @param polygon Polygon on which to modify the bathymetry
        """
        if polygon is None:
            self.set_array("MODEL.BOTTOMELEVATION", bathy)
        else:
            coordx = self.get_array("MODEL.X")
            coordy = self.get_array("MODEL.Y")
            for i, value in enumerate(bathy):
                if is_in_polygon(coordx[i], coordy[i], polygon):
                    self.set("MODEL.BOTTOMELEVATION",
                             value, i=i)

        return

    def __del__(self):
        """
        Destructor
        """
        Telemac2d._instanciated = False
