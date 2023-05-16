"""
Contains the class TelemacFile
"""
from os import path
import math
import re
import logging
from datetime import datetime

import numpy as np
from scipy.spatial import cKDTree
from scipy import interpolate
import matplotlib.tri as mtri

from utils.exceptions import TelemacException
from utils.polygon import points_in_poly
from telapy.api.hermes import HermesFile, ndp2elem, elem2str
from data_manip.extraction.linspace import \
        linspace_poly, curvilinear_abscissa

SPE_NAME = re.compile(r'F(?P<spe>.[0-9]{1,4})[ ]?PT2D(?P<point>[0-9]{6})')

class TelemacFile(HermesFile):
    """
    Class to extract data from a TelemacFile
    """

    def __init__(self, file_name, bnd_file=None, log_lvl='INFO', fformat=None,
                 access='r'):
        """
        Initialisation of a file reader

        @param file_name (str) Name of the mesh file
        @param bnd_file (str) Name of the boundary file (default None)
        @param log_lvl (str) Level of log information
        @param fformat (str) Format of the file if not given .med -> MED
        otherwise SERAFIN
        @param access (str) 'r' for read 'w' for write 'rw' for readwrite
        """

        # Identifying format from file extension if fformat not given
        if fformat is None:
            _, ext = path.splitext(file_name)
            if ext == '.med':
                fformat = 'MED'
            else:
                fformat = 'SERAFIN'

        HermesFile.__init__(self, file_name,
                            boundary_file=bnd_file,
                            access=access,
                            fformat=fformat,
                            log_lvl=log_lvl)

        self.interp_method = "matplotlib.LinearTri"
        self._title = None
        self._datetime = None

        self._ndim = None
        self._nelem3 = None
        self._npoin3 = None
        self._ndp3 = None
        self._nplan = None
        self._nelem2 = None
        self._npoin2 = None
        self._ndp2 = None
        self._meshx = None
        self._meshy = None
        self._meshz = None
        self._x_orig = None
        self._y_orig = None
        self._ikle3 = None
        self._ikle2 = None

        self._ipob3 = None
        self._ipob2 = None
        self._nptfr = None
        self._nelebd = None
        self._ikle_bnd = None
        self._nbor = None
        self._bnd_info = None

        self._nptir = None
        self._knolg = None

        self._ntimestep = None
        self._nvar = None
        self._varnames = None
        self._varunits = None
        self._values = None
        self._times = None

        self._tri = None
        self.tree = None
        self.neighbours = None
        self.edges = None

    def __del__(self):
        # To handle the case were we crashed in the hermes initialisation
        if "_title" in self.__dict__:
            del self._title
            del self._datetime

            del self._ndim
            del self._nelem3
            del self._npoin3
            del self._ndp3
            del self._nplan
            del self._nelem2
            del self._npoin2
            del self._ndp2
            del self._meshx
            del self._meshy
            del self._meshz
            del self._ikle3
            del self._ikle2

            del self._ipob3
            del self._ipob2
            del self._nptfr
            del self._nelebd
            del self._ikle_bnd
            del self._nbor
            del self._bnd_info

            del self._nptir
            del self._knolg

            del self._ntimestep
            del self._nvar
            del self._varnames
            del self._varunits
            del self._values
            del self._times

            del self._tri
            del self.tree
            del self.neighbours
            del self.edges

    @property
    def title(self):
        """
        Returns title value
        """
        if self._title is None:
            self._title = self.get_mesh_title()

        return self._title

    @title.setter
    def title(self, title):
        """ Setting title value """
        self._title = title

    @property
    def nvar(self):
        """
        Returns nvar value
        """
        if self._nvar is None:
            self._nvar = self.get_data_nvar()

        return self._nvar

    @nvar.setter
    def nvar(self, value):
        """ Setting nvar value """
        self._nvar = value

    @property
    def varnames(self):
        """
        Returns varnames value
        """
        if self._varnames is None:
            tmp_varnames, tmp_varunits = self.get_data_var_list()
            # Removing spaces at the end of the names/units
            self._varnames = [var.strip() for var in tmp_varnames]
            self._varunits = [var.strip() for var in tmp_varunits]

            del tmp_varnames
            del tmp_varunits

        return self._varnames

    @varnames.setter
    def varnames(self, value):
        """ Setting varnames value """
        self._varnames = value

    @property
    def varunits(self):
        """
        Returns varnames value
        """
        if self._varunits is None:
            self._varnames, self._varunits = self.get_data_var_list()

        return self._varunits

    @varunits.setter
    def varunits(self, value):
        """ Setting varunits value """
        self._varunits = value

    @property
    def datetime(self):
        """
        Returns datetime value
        """
        if self._datetime is None:
            self._datetime = self.get_mesh_date()

        return self._datetime

    @datetime.setter
    def datetime(self, value):
        """ Setting datetime value """
        self._datetime = value

    @property
    def ndim(self):
        """
        Return the dimnesion of the mesh
        """
        if self._ndim is None:
            self._ndim = self.get_mesh_dimension()

        return self._ndim

    @property
    def nelem3(self):
        """
        Returns number of 3d elements
        """
        if self._nelem3 is None:
            self._nelem3 = self.get_mesh_nelem()

        return self._nelem3

    @property
    def npoin3(self):
        """
        Returns number of 3d points
        """
        if self._npoin3 is None:
            self._npoin3 = self.get_mesh_npoin()

        return self._npoin3

    @property
    def ndp3(self):
        """
        Returns number of points in a 3d element
        """
        if self._ndp3 is None:
            self._ndp3 = self.get_mesh_npoin_per_element()

        return self._ndp3

    @property
    def nplan(self):
        """
        Returns the number of horizontal planes
        """
        if self._nplan is None:
            self._nplan = max(1, self.get_mesh_nplan())

        return self._nplan

    @property
    def nelem2(self):
        """
        Returns the number of 2d elements
        """
        if self._nelem2 is None:
            if self.nplan > 1:
                self._nelem2 = self.nelem3//(self._nplan-1)
            else:
                self._nelem2 = self.nelem3

        return self._nelem2

    @property
    def npoin2(self):
        """
        Returns the number of 2d points
        """
        if self._npoin2 is None:
            if self.nplan > 1:
                self._npoin2 = self.npoin3//self.nplan
            else:
                self._npoin2 = self.npoin3

        return self._npoin2

    @property
    def ndp2(self):
        """
        Returns the number of points per 2d element
        """
        if self._ndp2 is None:
            if self.nplan > 1:
                self._ndp2 = self.ndp3 // 2
            else:
                self._ndp2 = self.ndp3

        return self._ndp2

    @property
    def ikle3(self):
        """
        Returns the connectivity table for 3d elements
        """
        if self._ikle3 is None:
            self._ikle3 = self.get_mesh_connectivity()

        return self._ikle3

    @property
    def ikle2(self):
        """
        Returns the connectivity table for 2d elements
        """
        if self._ikle2 is None:
            if self.nplan > 1:
                self._ikle2 = np.compress(np.repeat([True, False], self.ndp2),
                                          self.ikle3[0:self.nelem2], axis=1)
            else:
                self._ikle2 = self.ikle3

        return self._ikle2

    @property
    def nptfr(self):
        """
        Returns the number of boundary points
        """
        if self._nptfr is None:
            if self.boundary_file is not None:
                self._nptfr = self.get_bnd_npoin()
            else:
                self._nptfr = 0

        return self._nptfr

    @property
    def nelebd(self):
        """
        Returns the number of boundary elements
        """
        if self._nelebd is None:
            if self.boundary_file is not None:
                self._nelebd = self.get_bnd_npoin()
            else:
                raise TelemacException(
                    "Can not read nelebd no boundary file was given")

        return self._nelebd

    @property
    def ikle_bnd(self):
        """
        Returns the number of boundary elements
        """
        if self._ikle_bnd is None:
            if self.boundary_file is not None:
                self._ikle_bnd = self.get_bnd_connectivity()
            else:
                raise TelemacException(
                    "Can not read ikle_bnd no boundary file was given")

        return self._ikle_bnd

    @property
    def ipob3(self):
        """
        Returns the ipobo array for 3d points
        """
        if self._ipob3 is None:
            if self.boundary_file is not None:
                self._ipob3 = self.get_bnd_ipobo()
            else:
                self._ipob3 = np.zeros((self.npoin3))

        return self._ipob3

    @property
    def ipob2(self):
        """
        Returns the ipobo array for 2d points
        """
        if self._ipob2 is None:
            if self.nplan > 1:
                self._ipob2 = self.ipob3[0:self._npoin2]
            else:
                self._ipob2 = self.ipob3

        return self._ipob2

    @property
    def meshx(self):
        """
        Returns x coordinates value
        """
        if self._meshx is None:
            self._meshx = self.get_mesh_coord(1)

        return self._meshx

    @property
    def meshy(self):
        """
        Returns y coordinates value
        """
        if self._meshy is None:
            self._meshy = self.get_mesh_coord(2)

        return self._meshy

    @property
    def meshz(self):
        """
        Returns y coordinates value
        """
        if self._meshz is None:
            if self.ndim == 3:
                self._meshz = self.get_mesh_coord(3)
            else:
                self._meshz = None

        return self._meshy

    @property
    def x_orig(self):
        """
        Returns y coordinates value
        """
        if self._x_orig is None:
            self._x_orig, self._y_orig = self.get_mesh_orig()

        return self._x_orig

    @property
    def y_orig(self):
        """
        Returns y coordinates value
        """
        if self._y_orig is None:
            self._x_orig, self._y_orig = self.get_mesh_orig()

        return self._y_orig

    @property
    def nptir(self):
        """
        Returns the number of interface points
        """
        if self._nptir is None:
            self._nptir = self.get_mesh_nptir()

        return self._nptir

    @property
    def knolg(self):
        """
        Return the local to global numbering array
        """
        if self._knolg is None:
            if self.nptir > 0:
                self._knolg = self.get_mesh_l2g_numbering()
            else:
                self._knolg = np.zeros((self.npoin3), dtype=np.int32)

        return self._knolg

    @property
    def ntimestep(self):
        """
        Returns the number of records in the file
        """
        if self._ntimestep is None:
            self._ntimestep = self.get_data_ntimestep()

        return self._ntimestep

    @property
    def tri(self):
        """
        Returns matplotlib triangulation of the 2d elements
        """
        if self._tri is None:
            self._tri = mtri.Triangulation(self.meshx[:self.npoin2],
                                           self.meshy[:self.npoin2],
                                           self.ikle2)

        return self._tri

    @property
    def times(self):
        """
        Returns a list of the times in the file
        """
        if self._times is None:
            self._times = np.zeros((self.ntimestep), dtype=np.float64)
            for record in range(self.ntimestep):
                self._times[record] = self.get_data_time(record)

        return self._times

    def set_kd_tree(self, reset=False):
        """
        Builds a KDTree (improves search of neighbours)

        @param reset (boolean) Force reset of tree
        """
        if reset or self.tree is None:
            isoxy = np.column_stack((np.sum(self.meshx[self.ikle2],
                                            axis=1)/3.0,
                                     np.sum(self.meshy[self.ikle2],
                                            axis=1)/3.0))
            self.tree = cKDTree(isoxy)

    def set_mpl_tri(self, reset=False):
        """
        Build neighbours from matplotlib

        @param reset (boolean) Force computing neighbours
        """
        if reset or self.neighbours is None or self.edges is None:
            # from matplotlib.tri import Triangulation
            mpltri = self.tri.get_cpp_triangulation()
            self.neighbours = mpltri.get_neighbors()
            self.edges = mpltri.get_edges()


    def get_z_name(self):
        """
        Return name of the variable containing the Z elevation
        If not found return None

        @returns (str) name of the variable containing the z elevation
        """
        if "COTE Z" in self.varnames:
            return "COTE Z"
        if "ELEVATION Z" in self.varnames:
            return "ELEVATION Z"

        return None


    #############################################
    #
    # Tools
    #
    #############################################

    def get_closest_record(self, time):
        """
        Get the record closest to a given time

        @param time (float) Time for which we seek the record

        @returns (int) The record
        """
        dist = 10000000
        record = -1
        for i, itime in enumerate(self.times):
            if (abs(time-itime)) < dist:
                dist = abs(time-itime)
                record = i

        return record

    def get_closest_node(self, point, plane=None):
        """
        Return the closest node to a given point If the mesh is 3d it will
        search for record 0

        @param point (np.array) coordinates of the point
        @param plane (int) (Only for a 3d mesh) If given will look for the
        closest node on the given plane point should be [x, y]

        @returns (int) Number of the node
        """
        node = -1
        best_dist = 1.e100

        if len(point) == 3:
            # Seaching in 3d mesh
            meshz = self.get_data_value(self.get_z_name(), 0)
            for i in range(self.npoin3):
                dist = (self.meshx[i] - point[0])**2 + \
                       (self.meshy[i] - point[1])**2 + \
                       (meshz[i] - point[2])**2

                if dist < best_dist:
                    best_dist = dist
                    node = i

        elif len(point) == 2:
            if plane is None:
                # Searching in a 2d mesh
                for i in range(self.npoin2):
                    dist = (self.meshx[i] - point[0])**2 + \
                           (self.meshy[i] - point[1])**2

                    if dist < best_dist:
                        best_dist = dist
                        node = i
            else:
                # Searching in a given plane for the closest node
                for i in range(plane*self.npoin2, (plane+1)*self.npoin2):
                    dist = (self.meshx[i] - point[0])**2 + \
                           (self.meshy[i] - point[1])**2

                    if dist < best_dist:
                        best_dist = dist
                        node = i

        else:
            raise TelemacException(
                "Point should be 2d or 3d: {}".format(point))

        return node

    def discretize_polyline(self, polyline):
        """
        Will return the number of point for each polyline segment taking the
        minimum mesh resolution as step

        @param polyline (list) List of points of the polyline

        @returns (list) List of discretisation for each segment
        """
        discret = []
        # ~~> Calculate the minimum mesh resolution
        dxy = math.sqrt(min(np.square(np.sum(np.fabs(
            self.meshx[self.ikle2] -
            self.meshx[np.roll(self.ikle2, 1)]),
                                             axis=1)/3.0) +
                            np.square(np.sum(np.fabs(
                                self.meshy[self.ikle2] -
                                self.meshy[np.roll(self.ikle2, 1)]),
                                             axis=1)/3.0)))
        for i in range(len(polyline)-1):
            dio = math.sqrt(sum(np.square(np.array(polyline[i])
                                          - np.array(polyline[i+1]))))
            discret.append(int(dio/dxy))

            del dio

        return discret

    #############################################
    #
    # Computing boundary information
    #
    #############################################
    @property
    def nbor(self):
        """
        Returns the boundary numbering
        """
        if self._nbor is None:
            if self.boundary_file is not None:
                self._nbor = self.get_bnd_numbering()
            else:
                raise TelemacException(
                    "Can not read nbor no boundary file was given")

        return self._nbor

    @property
    def bnd_info(self):
        """
        Get information for the boudnary file
        """
        if self._bnd_info is None:
            if self.boundary_file is not None:
                self._bnd_info = self.get_bnd_value()
            else:
                raise TelemacException(
                    "Can not read bnd_info no boundary file was given")

        return self._bnd_info

    def get_bnd_info(self):
        """
        Get boundary condition type of nodes
        """
        nbor = self.nbor
        lihbor, liubor, livbor, _, _, _, _, \
            litbor, _, _, _, _ = self.bnd_info

        return (nbor, lihbor, liubor, livbor, litbor)

    def get_liq_bnd_info(self):
        """
        Returns info on the liquid boundaries

        @returns
        """
        try:
            import _api as api
        except ImportError as xcpt:
            raise TelemacException(
                "Could not load the telemac api.\n"
                "They are mandatory for this function\n"+str(xcpt))

        ikles = self.ikle2.reshape(self.nelem2*3) + 1
        ndim = self.get_mesh_dimension()
        identify_liq_bnd = api.api_interface.identify_liq_bnd
        nbor, liubor, lihbor, _, _ = self.get_bnd_info()
        # Switching to fortran numbering
        nbor += 1

        coords = np.zeros((2, self.npoin2), dtype=np.float)
        coords[0, :] = self.meshx
        coords[1, :] = self.meshy
        coord = coords.reshape(self.npoin2*2)

        # Not using nelbor, ifabor, kp1bor (for now)
        _, _, _, numliq = identify_liq_bnd(
            ikles, ndim, liubor, lihbor, nbor,
            coord, False)

        nbor -= 1

        return nbor, numliq
    #############################################
    #
    # Interpolation function
    #
    #############################################
    def interpolate_mesh_on_points(self, points, values):
        """

        Interpolate mesh on points given value

        @param points (list) List og points for which to interpolate
        @param values (np.darray) values for each mesh points
        @param method (string) Interpolation methof to use

        """
        if self.interp_method == "scipy.LinearND":
            grid_pt = np.column_stack((self.meshx[:self.npoin2],
                                       self.meshy[:self.npoin2]))
            data_interp = interpolate.LinearNDInterpolator(grid_pt, values)
            pts = np.asarray([(pt[0], pt[1]) for pt in points])
            res = data_interp(pts)
        elif self.interp_method == "matplotlib.LinearTri":
            data_interp = mtri.LinearTriInterpolator(self.tri, values)
            pt_x = [pt[0] for pt in points]
            pt_y = [pt[1] for pt in points]
            res = data_interp(pt_x, pt_y)
        else:
            raise TelemacException(\
               "Unknown interpolation method: {}\n"
               "Avaialable methods: matplotlib.LinearTri, scipy.LinearND"\
               .format(self.interp_method))

        return res

    def interpolate_xyz_on_mesh(self, xyz, kind='nearest', poly=None,
                                loc='inside', init=None,
                                px=None, py=None, pe=0.25,\
                                **kwargs):
        """
        Interpolate xyz information on the TelemacFile mesh
        Additional position argument will be passed to the interpolator Interp2D

        @param
        """
        from pretel.interpolator import Interp2D, Interp2DChunk
        # Mesh coordonate on which to interpolate data
        if poly is not None:
            coord_x = self.meshx
            coord_y = self.meshy
            points = np.column_stack((coord_x, coord_y))
            inside = points_in_poly(points, poly)
            if loc == "inside":
                inloc = inside
            elif loc == "outside":
                inloc = np.logical_not(inside)
            coord_x = coord_x[inloc]
            coord_y = coord_y[inloc]
        else:
            coord_x = self.meshx
            coord_y = self.meshy

        # Interpolator
        if px is not None and py is not None:
            interp = Interp2DChunk(xyz, kind=kind, \
                                   px=px, py=py, pe=pe, **kwargs)
        else:
            interp = Interp2D(xyz, kind=kind, **kwargs)

        data = interp(coord_x, coord_y)

        # Interpolate data on mesh
        if poly is not None:
            if init is None:
                data_full = np.nan((self.npoin2))
            elif isinstance(init, (int, float)):
                data_full = np.ones((self.npoin2), dtype=np.float64)*init
            elif isinstance(init, (np.ndarray, list)):
                if len(init) != self.npoin2:
                    raise TelemacException(
                        "Init array should be size {} but is {}"
                        .format(self.npoin2, len(init)))
                data_full = init.copy()
            else:
                raise TelemacException(
                    "type {} is not handled for init".format(type(init)))
            data_full[inloc] = data
        else:
            data_full = data

        return data_full

    #############################################
    #
    # data extractions functions
    #
    #############################################
    #
    # Extractrion for a given record
    #
    def get_data_on_points(self, varname, record, points):
        """
        Extract values on points in telemac result file (2D or 3D)
        for the given variable for one record

        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of desired record to extract
        @param points list of numpy.array containing points of extraction

        @returns (numpy.array)
        """
        if len(np.shape(np.array(points))) != 2:
            raise TelemacException('Warning problem with the list of '
                                   'extraction points')
        # dimension of the computation result
        dim = np.shape(np.array(points))[1]
        if dim == 2:
            res = self._get_data_on_2d_points(varname, record, points)
        elif dim == 3:
            res = self._get_data_on_3d_points(varname, record, points)
        else:
            raise TelemacException('Warning problem with the dimension of '
                                   'extraction points')
        return res

    def _get_data_on_2d_points(self, varname, record, points):
        """
        Extract values on points in telemac-2d result file
        for the given variable for one record

        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of the desired record to extract
        @param points (list) list of points to extract

        @returns (numpy.array)
        """
        if self.get_mesh_dimension() != 2:
            raise TelemacException("Action possible only on 2d mesh")

        values = self.get_data_value(varname, record)
        if len(values) > self.npoin2:
            raise TelemacException('Warning the dimension of the result '
                                   'file is greater than 2')

        res = self.interpolate_mesh_on_points(points, values)

        return res

    def _get_data_on_3d_points(self, varname, record, points):
        """
        Extract values on points in telemac-3d result file
        for the given variable for one record

        @param points list of numpy.array containing points of
               extraction (x,y,z)
        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of desired record to extract

        @returns (numpy.array)
        """
        if self.get_mesh_dimension() != 3:
            raise TelemacException("Action possible only on 3d mesh")

        res = float('nan')*np.ones((len(points)), dtype=np.float64)
        for i, point in enumerate(points):
            elev = self.get_data_on_vertical_segment(
                self.get_z_name(), record, point[:-1])
            values = self.get_data_on_vertical_segment(
                varname, record, point[:-1])
            for plan in range(self.nplan-1):
                if elev[plan] <= point[-1] and point[-1] <= elev[plan+1]:
                    shz = (point[-1]-elev[plan])/max((elev[plan+1]
                                                      - elev[plan]), 1.e-6)
                    res[i] = (1.0-shz)*values[plan]+shz*values[plan+1]
        return res

    def get_data_on_polyline(self, varname, record, polyline_points,
                             discretized_number=None):
        """
        Extract values of points over time for the given variable
        for record

        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of the desired record to extract
        @param polyline_points (list) List of points defining the polyline
        @param discretized_number (list) list of number of discretized points
        on each polyline segment if None given will
            use self.discretize_polyline

        @returns (numpy.array, numpy.array, numpy.array)

        """
        if self.get_mesh_dimension() != 2:
            raise TelemacException("Action possible only on 2d mesh")

        if len(np.shape(np.array(polyline_points))) != 2:
            raise TelemacException('Warning problem with the list of '
                                   'extraction points')

        if discretized_number is None:
            discretized_number = self.discretize_polyline(polyline_points)

        # dimension of the computation result
        dim = np.shape(np.array(polyline_points))[1]
        if dim == 2:
            polygone_discretized_points = linspace_poly(polyline_points,
                                                        discretized_number)
            values_polylines = self.get_data_on_points(
                varname,
                record,
                polygone_discretized_points)
            abs_curv = curvilinear_abscissa(polygone_discretized_points)
        else:
            raise TelemacException('Warning the extraction on a polyline'
                                   ' is valid only in 2d')
        return polygone_discretized_points, abs_curv, values_polylines

    def get_data_on_horizontal_plane(self, varname, record, plane_number):
        """
        Extract values of one plane in telemac-3d result file
        for the given variable

        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of desired record
        @param plane_number (int) Number of desired plane

        @returns (numpy.array)

        """
        if self.get_mesh_dimension() != 3:
            raise TelemacException("Action possible only on 3d mesh")

        values = self.get_data_value(varname, record)
        if plane_number < 0 or plane_number >= self.nplan:
            raise TelemacException(
                'Wrong plane number {} should be in [0, {}]'
                .format(plane_number, self.nplan-1))
        start = plane_number*self.npoin2
        end = (plane_number+1)*self.npoin2
        extracted_values = values[start:end]

        return extracted_values

    def get_data_on_horizontal_slice(
            self, varname, record, zslices, nplanref=None):
        """
        Extract values of plan in telemac-3d result file for the given variable

        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of desired record
        @param zslices (numpy.array/list/int) Elevation of the slice or
        @param nplanref (int) Number of reference plane

        @returns (numpy.array)
        """
        if self.get_mesh_dimension() != 3:
            raise TelemacException("Action possible only on 3d mesh")

        if isinstance(zslices, (list, np.ndarray)):
            if isinstance(zslices, np.ndarray):
                if zslices.ndim > 1:
                    raise TelemacException('Warning the slice coordinate'
                                           'must be 1d')
            res = np.zeros(((self.npoin2), len(zslices)), dtype=np.float64)
            zslices_list = zslices
        elif isinstance(zslices, int):
            res = np.zeros(((self.npoin2), zslices), dtype=np.float64)
            zslices_list = [zslices]
        else:
            raise TelemacException('Unknown zslices type')

        zref = np.zeros((self.npoin2), dtype=np.float64)

        z_name = self.get_z_name()

        if z_name is None:
            raise TelemacException('Variable for elevation is missing')

        if nplanref is not None:
            zref = self.get_data_on_horizontal_plane(
                z_name, record, nplanref)
        values_elevation = self.get_data_value(z_name, record)
        values_elevation = values_elevation.reshape(self.nplan,
                                                    self.npoin2)
        values_var = self.get_data_value(varname, record)
        values_var = values_var.reshape(self.nplan, self.npoin2)

        for izs, zslice in enumerate(zslices_list):
            zslice = zref + zslice
            for j in range(self.npoin2):
                res[j, izs] = float('nan')
                for i in range(self.nplan-1):
                    if values_elevation[i, j] <= zslice[j] and \
                       zslice[j] <= values_elevation[i+1, j]:
                        shz = (zslice[i]-values_elevation[i, j]) /\
                              max((values_elevation[i+1, j]
                                   - values_elevation[i, j]), 1.0e-6)
                        res[j, izs] = (1.0-shz)*values_var[i, j]+shz *\
                            values_var[i+1, j]
                        break

        if isinstance(zslices, (list, np.ndarray)):
            return res
        if isinstance(zslices, int):
            return res[:, 0]

        return None

    def get_data_on_vertical_plane(self, varname, record, polyline_points,
                                   discretized_number=None):
        """
        Extract values of plan in telemac-3d result file for the given variable

        @param varname (string) Name of variable for which to extract data
        @param record (int) Number of desired record
        @param polyline_points (list) List of points defining the polyline
        @param discretized_number (list) List of number of discretized points
        on each polyline segment

        @returns (numpy.array)

        """
        if self.get_mesh_dimension() != 3:
            raise TelemacException("Action possible only on 3d mesh")

        if len(np.shape(np.array(polyline_points))) != 2:
            raise TelemacException('Warning problem with the list of '
                                   'extraction points')

        if discretized_number is None:
            discretized_number = self.discretize_polyline(polyline_points)

        dim = np.shape(np.array(polyline_points))[1]
        if dim == 2:
            nplan = self.nplan
            polygone_discretized_points = linspace_poly(
                polyline_points, discretized_number)
            npoly = len(polygone_discretized_points)
            values_polylines = np.zeros((npoly, nplan), dtype=np.float64)
            abs_curv = curvilinear_abscissa(polygone_discretized_points)
            for plan in range(self.nplan):
                values = self.get_data_on_horizontal_plane(
                    varname, record, plan)
                values_polylines[:, plan] = \
                        self.interpolate_mesh_on_points(\
                            polygone_discretized_points, values)
                del values
        else:
            raise TelemacException('Warning the extraction on a polyline'
                                   ' of 2d points')
        return polygone_discretized_points, abs_curv, values_polylines

    def get_data_on_vertical_segment(self, var_name, record, point):
        """
        Extract values for each plane of a 2d points in telemac-3d result file
        for the given variable

        @param var_name (string) Name of variable for which to extract data
        @param record (int) Number of desired record
        @param point (numpy.array) Point of extraction

        @returns (numpy.array)

        """
        if self.get_mesh_dimension() != 3:
            raise TelemacException("Action possible only on 3d mesh")

        if len(point) != 2:
            raise TelemacException('Warning the extraction point '
                                   'must be 2d')
        nplan = self.nplan
        res = np.zeros(nplan)
        for plan in range(self.nplan):
            values = self.get_data_on_horizontal_plane(
                var_name, record, plan)
            res[plan] = self.interpolate_mesh_on_points([point], values)
        return res

    #
    # Extractrion of a timeserie (extraction of data over all records)
    #
    def get_timeseries_on_nodes(self, varname, nodes):
        """
        Extract values of nodes over time for the given variable

        @param varname (string) Name of variable for which to extract data
        @param nodes (list) list of nodes to extract

        @returns (numpy.array) shape (len(nodes), self.ntimestep)

        """
        res = np.zeros((len(nodes), self.ntimestep), dtype=np.float64)
        for record in range(self.ntimestep):
            values = self.get_data_value(varname, record)
            res[range(len(nodes)), record] = values[nodes]

        return res

    def get_timeseries_on_points(self, varname, points):
        """
        Extract values of points over time for the given variable

        @param varname (string) Name of variable for which to extract data
        @param points (list) List of points to extract

        @returns (numpy.array)

        """
        res = np.zeros((len(points), self.ntimestep), dtype=np.float64)
        for record in range(self.ntimestep):
            res[:, record] = self.get_data_on_points(varname, record, points)
        return res

    def get_timeseries_on_polyline(self, varname, polyline_points,
                                   discretized_number=None):
        """
        Extract values of points over time for the given variable

        @param varname (string) Name of variable for which to extract data
        @param polyline_points (list) List of points defining the polyline
        @param discretized_number (list) List of number of discretized points
        on each polyline segment

        @returns (numpy.array, numpy.array, numpy.array) polygone discretised
        polygone_discretized_points, abs_curv, values_polylines

        """
        if self.get_mesh_dimension() != 2:
            raise TelemacException("Action possible only on 2d mesh")

        if len(np.shape(np.array(polyline_points))) != 2:
            raise TelemacException('Warning problem with the list of '
                                   'extraction points')

        if discretized_number is None:
            discretized_number = self.discretize_polyline(polyline_points)

        # dimension of the computation result
        dim = np.shape(np.array(polyline_points))[1]
        if dim == 2:
            polygone_discretized_points = linspace_poly(polyline_points,
                                                        discretized_number)
            values_polylines = self.get_timeseries_on_points(
                varname, polygone_discretized_points)
            abs_curv = curvilinear_abscissa(polygone_discretized_points)
        else:
            raise TelemacException('Warning the extraction on a polyline'
                                   ' is valid only in 2d')
        return polygone_discretized_points, abs_curv, values_polylines

    def get_timeseries_on_vertical_segment(self, varname, point):
        """
        Extract values of plan in telemac-3d result file for the given variable

        @param varname (string) Name of variable for which to extract data
        @param point (numpy.array) Point of extraction

        @returns (numpy.array)

        """
        if self.get_mesh_dimension() != 3:
            raise TelemacException("Action possible only on 3d mesh")

        res = np.zeros((self.nplan, self.ntimestep), dtype=np.float64)
        for record in range(self.ntimestep):
            res[:, record] = self.get_data_on_vertical_segment(
                varname, record, point)
        return res

    #############################################
    #
    # Spectrum specific functions
    #
    #############################################
    def is_a_spectrum_file(self):
        """
        Checking if the file is a spectrum file
        Criteria:
        - Quadrangles
        - All variables start with FPTS
        """

        is_spectrum = self.ndp3 == 4
        #
        is_spectrum = is_spectrum and \
            all([SPE_NAME.match(var) is not None for var in self.varnames])

        return is_spectrum

    def get_spectrum_freq(self):
        """
        Compute the list of frequencies
        This only works if the file is a tomawac spectrum file

        @returns (numpy.array, numpy.array) List of frequencies, List of
        frequencie steps
        """
        if not self.is_a_spectrum_file():
            raise TelemacException(
                "This file does not seem to be a spectrum file")

        nfreq = 0
        eps = 1e-6
        f_1 = 10e10
        f_2 = 10e10
        raisf = 0.
        for x, y in zip(self.meshx, self.meshy):
            if abs(x) <= eps and y >= 0.:
                nfreq += 1
                f_temp = y
                if f_temp < f_1:
                    f_2 = f_1
                    f_1 = f_temp
                elif f_temp < f_2:
                    f_2 = f_temp

        raisf = f_2/f_1

        freqs = [f_1 * raisf**i for i in range(nfreq)]

        dfreqs = np.zeros(nfreq, dtype=np.float64)

        auxi = (raisf - 1.)/2.
        dfreqs[0] = auxi*freqs[0]
        for i in range(1, nfreq-1):
            dfreqs[i] = auxi*(freqs[i] + freqs[i-1])

        dfreqs[-1] = auxi*freqs[-2]

        return np.array(freqs), dfreqs

    def get_list_spectrum_points(self):
        """
        Returns the list of spectrum points in the file
        """
        if not self.is_a_spectrum_file():
            raise TelemacException(
                "This file does not seem to be a spectrum file")

        points = []
        for var in self.varnames:
            proc = SPE_NAME.match(var)
            number = proc.group('point').lstrip('0')
            points.append(int(number))

        return points

    def get_spectrum_varname(self, point):
        """
        Return the variable associated to the spectrum point 'point'

        @param point (int) Point number

        @returns (string) Name of the variable
        """
        if not self.is_a_spectrum_file():
            raise TelemacException(
                "This file does not seem to be a spectrum file")

        spectrum_var = None
        # Getting the variable for point point
        for var in self.varnames:
            if "{:06d}".format(point) in var:
                spectrum_var = var
                break

        if spectrum_var is None:
            raise TelemacException("Could not find point {} in your variables:\
            \n{}".format(point, self.varnames))

        return spectrum_var

    def get_angular_dispersion(self, point, record, radian=False):
        """
        Return value of the angular dispersion

        @param point (int) number of the point for which
            we extract the spectrum
        @param record (int) Time record for which to extract
        @param radian (boolean) If true theta is built in radian otherwise in
        degree

        @returns (numpy.array, numpy.array) The frequencie list, The angular
        dispersion values
        """
        if not self.is_a_spectrum_file():
            raise TelemacException(
                "This file does not seem to be a spectrum file")

        spectrum_var = self.get_spectrum_varname(point)

        # Getting list of frequencies
        freqs, dfreqs = self.get_spectrum_freq()

        nfreq = len(freqs)
        ntheta = self.npoin2//nfreq

        # Reshaping to match nfreq*ntheta
        data = self.get_data_value(spectrum_var, record)\
                   .reshape((nfreq, ntheta))

        ang_disp = np.zeros(ntheta, dtype=np.float64)
        # Integration over frequencies
        for itheta in range(ntheta):
            for ifreq in range(nfreq):
                ang_disp[itheta] += data[ifreq, itheta]*dfreqs[ifreq]

        # Defining if we are in radian or degree
        if radian:
            val = 2*np.pi
        else:
            val = 360.

        # Building angles array
        theta = [i*val/ntheta for i in range(ntheta)]

        return theta, ang_disp

    def get_spectrum(self, point, record):
        """
        Return value of spectrum for a given point and record

        @param point (int) number of the point for which
            we extract the spectrum
        @param record (int) Time record for which to extract

        @returns (numpy.array, numpy.array) The frequencie list, The spectrum
        values
        """
        if not self.is_a_spectrum_file():
            raise TelemacException(
                "This file does not seem to be a spectrum file")

        spectrum_var = self.get_spectrum_varname(point)

        # Getting list of frequencies
        freqs, _ = self.get_spectrum_freq()

        nfreq = len(freqs)
        ntheta = self.npoin2//nfreq

        # Reshaping to match nfreq*ntheta
        data = self.get_data_value(spectrum_var, record)\
                   .reshape((nfreq, ntheta))

        # Integration over angles
        spectrum = np.sum(data, axis=1) * 2*np.pi/ntheta

        return freqs, spectrum
    ###
    # Writing functions
    ###
    def add_variable(self, varname, varunit):
        """
        Adding new variable to file (If no time step in file will create one
        with time 0.0) in the class

        @param varname (str) Name of the variable
        @param varname (str) Unit of the variable
        @param values (np.array) If given values for that variable that will be
        set on all time step if not given variable will be set to zero
        """
        if self._varnames is None:
            self._varnames = []
            self._varunits = []
            self._nvar = 0

        self._varnames.append(varname)
        self._varunits.append(varunit)
        self._nvar += 1

        if self._ntimestep is None:
            self._ntimestep = 1
            self._times = np.array([0.0], dtype=np.float64)

        if self._values is None:
            self._values = np.zeros((self._ntimestep, self._nvar, self._npoin3),
                                    dtype=np.float64)
        else:
            data = np.zeros((self._ntimestep, 1, self._npoin3))
            self._values = np.append(self._values, data, axis=1)

    def add_time_step(self, time):
        """
        Adding a new time step if some variables are set expand values array in
        the class

        @param time (float)
        """
        if self._ntimestep is None:
            self._ntimestep = 1
            self._times = np.array([time], dtype=np.float64)
        else:
            self._times = np.append(self.times, [time], axis=0)
            self._ntimestep += 1

        if self._values is not None:
            data = np.zeros((1, self._nvar, self._npoin3))
            self._values = np.append(self._values, data, axis=0)

    def add_header(self, title, date=None):
        """
        Adding header data in the class

        @param title (str) Title of the file
        @date date (np.array shape (6)) Date and time
        """
        self._title = title
        self._datetime = date

    def add_mesh(self, x, y, ikle, z=None, nplan=0, orig=(0, 0)):
        """
        Adding mesh information in the class

        @param x (list/np.array) X coordinates of the mesh
        @param y (list/np.array) Y coordinates of the mesh
        @param z (list/np.array) Z coordinates of the mesh (default)
        @param ikle (list/np.array (nelem, ndp)) Connectivity of mesh
        @param nplan (int) Number of plane (default 0)
        @param orig (in, int) Origin of the mesh (default (0, 0))
        """

        # Points info
        self._meshx = np.array(x, dtype=np.float64)
        self._meshy = np.array(y, dtype=np.float64)
        if z is not None:
            self._meshz = z
        self._npoin3 = self._meshx.shape[0]

        self._x_orig, self._y_orig = orig

        # Connectivity info
        self._ikle3 = np.array(ikle, dtype=np.int32)
        self._nelem3, self._ndp3 = self._ikle3.shape

        self._nplan = nplan

        self._ndim = 2 if nplan <= 1 else 3

        self.typ_elem = ndp2elem(self._ndp3, self._ndim)

        # variable initialised to 0 so that the class is valid for a write
        self._nptfr = 0
        self._nptir = 0
        self._ipob3 = np.zeros((self._npoin3))
        self._knolg = np.zeros((self._npoin3))

    def add_bnd(self, ikle, lihbor=None, liubor=None, livbor=None,
                hbor=None, ubor=None, vbor=None, chbord=None,
                litbor=None, tbor=None, atbor=None, btbor=None,
                color=None):
        """
        Adding boundary information in the class

        @param ikle (list/np.array) boudnary connectivity
        """

        self._ikle_bnd = np.array(ikle)

        if self._ikle_bnd.ndim == 1:
            ndp = 1
            self._nptfr = self._ikle_bnd.shape[0]
            self._nelebd = self._nptfr
            self._nbor = self._ikle_bnd.copy()
        else:
            self._nelebd, ndp = self._ikle_bnd.shape
            self._nptfr = self._nelebd -1
            self._nbor = self._ikle_bnd[:, 0]

        self.typ_bnd_elem = ndp2elem(ndp, self._ndim)

        self._ipob3[:] = 0
        self._ipob3[self._nbor] = self._nbor

        if lihbor is None:
            lihbor = np.ones((self._nelebd)) * 2
        if liubor is None:
            liubor = np.ones((self._nelebd)) * 2
        if livbor is None:
            livbor = np.ones((self._nelebd)) * 2
        if hbor is None:
            hbor = np.zeros((self._nelebd), dtype=np.float64)
        if ubor is None:
            ubor = np.zeros((self._nelebd), dtype=np.float64)
        if vbor is None:
            vbor = np.zeros((self._nelebd), dtype=np.float64)
        if chbord is None:
            chbord = np.zeros((self._nelebd), dtype=np.float64)
        if litbor is None:
            litbor = np.ones((self._nelebd)) * 2
        if tbor is None:
            tbor = np.zeros((self._nelebd), dtype=np.float64)
        if atbor is None:
            atbor = np.zeros((self._nelebd), dtype=np.float64)
        if btbor is None:
            btbor = np.zeros((self._nelebd), dtype=np.float64)
        if color is None:
            color = np.array([i+1 for i in range(self._nelebd)])

        self._bnd_info = lihbor, liubor, livbor, hbor, ubor, vbor, chbord, \
            litbor, tbor, atbor, btbor, color


    def add_data_value(self, varname, record, values):
        """
        Adding values for a given variable and record in the class
        """
        var_idx = self.varnames.index(varname)


        # Define data to write from values type
        if isinstance(values, (int, float)):
            data = np.ones((self._npoin3), dtype=np.float64) * values
        elif isinstance(values, (np.ndarray, list)):
            if len(values) != self._npoin3:
                raise TelemacException(
                    "Values is of list/numpy type it is {} long is should be {}"
                    .format(len(values), self._npoin3))
            data = np.array(values)
        elif callable(values):
            # Check the function as two arguments
            import inspect
            if len(inspect.getargspec(values)[0]) != 2:
                raise TelemacException(
                    "values function must take two arguments")
            data = values(self._meshx, self._meshy)
        else:
            raise TelemacException("type {} is not handled for values"
                                   .format(type(values)))

        self._values[record, var_idx, :] = data


    def set_data_on_polygon(self, varname, record, values,
                            poly, loc='inside'):
        """
        Set data value on triangular mesh within/without a polygon

        @param varname (str) Name of the variable to write
        @param record (int) record for which to write (if < 0 uses ntimespte+record)
        @param values (float/list/numpy/function) Values to set it can be:
                       - a constant
                       - a array of npoin3 size
                       - a function of x and y that will be applied to the coordiantes of the points
        """
        points = np.column_stack((self.tri.x, self.tri.y))

        # Identify points in polygon
        inside = points_in_poly(points, poly)

        if loc == "inside":
            inloc = inside
        elif loc == "outside":
            inloc = np.logical_not(inside)

        coord_x = self.meshx[inloc]
        coord_y = self.meshy[inloc]

        data_full = self.get_data_value(varname, record)
        # check that values is of proper size
        if isinstance(values, (int, float)):
            data = np.ones((coord_x.shape[0]), dtype=np.float64) * values
        elif callable(values):
            # Check the function as two arguments
            import inspect
            if len(inspect.getargspec(values)[0]) != 2:
                raise TelemacException(
                    "values function must take two arguments")
            data = values(coord_x, coord_y)
        else:
            raise TelemacException("type {} is not handled for values"
                                   .format(type(values)))

        data_full[inloc] = data

        self.set_data_value(varname, record, data_full)

    def set_data_value(self, varname, record, values):
        """
        Set values for a given variable and record

        @param varname (str) Name of the variable to write
        @param record (int) record for which to write (if < 0 uses
        ntimestee+record)
        @param values (float/list/numpy/function) Values to set it can be:
                       - a constant
                       - a array of npoin3 size
                       - a function of x and y that will be applied to the
                         coordiantes of the points
        """

        # Check that variable is in file
        if varname not in self.varnames:
            raise TelemacException(
                "{} is not the file variables:\n{}"
                .format(varname, self.varnames))
        var_idx = self.varnames.index(varname)

        # Check that record is in file
        if record < 0:
            record = self.ntimestep + record
        if record < 0 or record >= self.ntimestep:
            raise TelemacException(
                "Record {} is not within [0, {}]"
                .format(record, self.ntimestep))

        # Define data to write from values type
        if isinstance(values, (int, float)):
            data = np.ones((self.npoin3), dtype=np.float64) * values
        elif isinstance(values, (np.ndarray, list)):
            if len(values) != self.npoin3:
                raise TelemacException(
                    "Values is of list/numpy type it is {} long is should be {}"
                    .format(len(values), self.npoin3))
            data = np.array(values)
        elif callable(values):
            # Check the function as two arguments
            import inspect
            if len(inspect.getargspec(values)[0]) != 2:
                raise TelemacException(
                    "values function must take two arguments")
            data = values(self.meshx, self.meshy)
        else:
            raise TelemacException("type {} is not handled for values"
                                   .format(type(values)))

        self.add_data(varname, self.varunits[var_idx], record,
                      self.times[record], var_idx == 0, data)

    ###
    # Loading/Dumping functions
    ###

    def read(self, src=None):
        """
        This function will read all the information within a
        telemac file and store it in its variables

        @param src (TelemacFile) Object from which to read data if none given
        will use self
        """
        if src is None:
            src = self
        self.read_mesh(src)
        self.read_data(src)

    def read_mesh(self, src=None):
        """
        This function will read all the information within the
        telemac file and store it in its variables

        @param src (TelemacFile) Object from which to read data if none given
        will use self
        """
        if src is None:
            src = self
        self.logger.debug("Reading mesh information from file %s",
                          src.file_name)

        self._title = src.title
        self._datetime = src.datetime

        self._ndim = src.ndim

        # copying mesh quantities
        self._npoin3 = src.npoin3
        self._nelem3 = src.nelem3
        self._ndp3 = src.ndp3
        self._nplan = src.nplan
        self.typ_elem = src.typ_elem

        # Copying mesh coordiantes
        self._meshx = src.meshx
        self._meshy = src.meshy
        self._meshz = src.meshz
        self._x_orig = src.x_orig
        self._y_orig = src.y_orig

        # Copying connectivity
        self._ikle3 = src.ikle3

        # Parallel interface information
        self._nptir = src.nptir
        self._knolg = src.knolg

        # Boundary information
        # nptfr and ipob3 are read reagrdless of presence of boundary file
        # As they are need in serafin format
        self._nptfr = src.nptfr
        self._ipob3 = src.ipob3
        if self.boundary_file is not None:
            self.typ_bnd_elem = src.typ_bnd_elem
            self._nelebd = src.nelebd
            self._bnd_info = src.bnd_info
            self._ikle_bnd = src.ikle_bnd
            self._nbor = src.nbor

    def read_data(self, src=None):
        """
        Read information on fields (variables, records, value for each variable
        and record)
        Warning this function can be very memory consuming if you have lots of
        records and variables

        @param src (TelemacFile) Object from which to read data if none given
        will use self
        """
        if src is None:
            src = self
        self.logger.debug("Reading data information from file %s",
                          src.file_name)


        self._ntimestep = src.get_data_ntimestep()
        self._nvar = src.get_data_nvar()
        self._varnames, self._varunits = src.get_data_var_list()

        self._times = np.zeros((self._ntimestep), dtype=np.float64)
        self._values = np.zeros((self._ntimestep, self._nvar, self._npoin3),
                                dtype=np.float64)

        for i in range(self._ntimestep):
            self._times[i] = src.get_data_time(i)
            for j in range(self._nvar):
                self._values[i, j] = src.get_data_value(self._varnames[j], i)

    def write(self):
        """
        Writting data from class into file
        """
        # Header part
        self.logger.debug("Writting header information from class in file %s",
                          self.file_name)
        # Checking that variables are properly sets
        for variable in ['title', 'varnames', 'varunits', 'nvar']:
            if getattr(self, "_"+variable) is None:
                raise TelemacException("Missing {} in class".format(variable))

        # Checking dimensions of varnames and varunits
        if len(self._varnames) != self._nvar:
            raise TelemacException(
                "Error in varnames we have {} variables and {} names"
                "\n varnames: {}"
                .format(self._nvar, len(self._varnames), self._varnames))
        if len(self._varunits) != self._nvar:
            raise TelemacException(
                "Error in varnames we have {} variables and {} units"
                "\n varunits: {}"
                .format(self._nvar, len(self._varunits), self._varunits))

        self.set_header(self._title, self._nvar, self._varnames,
                        self._varunits)

        # Mesh part
        date2 = np.zeros((3), dtype=np.int32)
        time2 = np.zeros((3), dtype=np.int32)
        date2[0] = self._datetime[0]
        date2[1] = self._datetime[1]
        date2[2] = self._datetime[2]
        time2[0] = self._datetime[3]
        time2[1] = self._datetime[4]
        time2[2] = self._datetime[5]

        self.logger.debug("Writting mesh information from class in file %s",
                          self.file_name)
        # Checking that variables are properly sets
        for variable in ['ndim', 'ndp3', 'nptfr', 'nptir', 'nelem3', 'npoin3',
                         'ikle3', 'ipob3', 'knolg', 'meshx', 'meshy', 'nplan',
                         'datetime']:
            if getattr(self, "_"+variable) is None:
                raise TelemacException("Missing {} in class".format(variable))

        self.set_mesh(self._ndim, self.typ_elem, self._ndp3, self._nptfr,
                      self._nptir, self._nelem3, self._npoin3,
                      self._ikle3, self._ipob3, self._knolg,
                      self._meshx, self._meshy, self._nplan, date2,
                      time2, self._x_orig, self._y_orig, self._meshz)

        # Boundary part
        if self.boundary_file is not None:
            self.logger.debug("Writting bnd information from class in file %s",
                              self.file_name)
            # Checking that variables are properly sets
            for variable in ['nelebd', 'ikle_bnd', 'bnd_info', 'nbor']:
                if getattr(self, "_"+variable) is None:
                    raise TelemacException(
                        "Missing {} in class".format(variable))

            lihbor, liubor, livbor, hbor, ubor, vbor, chbord, \
                litbor, tbor, atbor, btbor, color = self._bnd_info

            self.set_bnd(self.typ_bnd_elem, self._nelebd, self._ikle_bnd,
                         lihbor, liubor, livbor, hbor, ubor, vbor, chbord,
                         litbor, tbor, atbor, btbor, color)

        # Data part
        self.logger.debug("Writting data information from class in file %s",
                          self.file_name)

        # Checking that variables are properly sets
        for variable in ['ntimestep', 'times', 'nvar', 'values']:
            if getattr(self, "_"+variable) is None:
                raise TelemacException("Missing {} in class".format(variable))

        # Chacking dimensions of values and times
        if self._values.shape != (self._ntimestep, self._nvar, self._npoin3):
            raise TelemacException(
                "Error in shape of values (ntimestep, nvar, npoin3):"
                "\nvalues is {} and should be {}"
                .format(self._values.shape,
                        (self._ntimestep, self._nvar, self._npoin3)))

        if self._times.shape != (self._ntimestep,):
            raise TelemacException(
                "Error in shape of times (ntimestep):"
                "\ntimes is {} and should be {}"
                .format(self._times.shape,
                        (self._ntimestep,)))

        for i in range(self._ntimestep):
            time = self._times[i]
            for j in range(self._nvar):
                self.add_data(self._varnames[j], self._varunits[j],
                              time, i, j == 0, self._values[i, j])

    def print_info(self, full=False):
        """
        Printing mesh information

        @param full (bool) If True display data informations
        """
        print('~> Generic info')
        print()

        # TODO: print full title or without format ?
        print("  - Title: {}".format(self.title[:-8].strip()))
        date = self.datetime
        if any(date != 0):
            date = datetime(date[0], date[1], date[2],
                            date[3], date[4], date[5])
            print("  - Date: {}".format(date))
        else:
            print("  - No date in file")

        print("  - Format: {}".format(self.fformat.decode('utf-8')))
        #TODO: Add endianess info

        print()
        print('~> Mesh info')
        print()

        print("  - Number of dimensions: {}".format(self.ndim))
        print("  - Element type: {}".format(elem2str(self.typ_elem)))

        print("  - Number of points: {}".format(self.npoin3))
        print("  - Number of elements: {}".format(self.nelem3))
        print("  - Number of points per element: {}".format(self.ndp3))
        print("  - Number of planes: {}".format(self.nplan))

        if self.ndim == 3:
            print("  - Number of points: {}".format(self.npoin2))
            print("  - Number of elements: {}".format(self.nelem2))
            print("  - Number of points per element: {}".format(self.ndp2))

        print()
        print("  +> Coordinates")
        print()
        print("    - X offset, Y offset = {}, {}"
              .format(self.x_orig, self.y_orig))
        x_min = np.min(self.meshx)
        y_min = np.min(self.meshy)
        x_max = np.max(self.meshx)
        y_max = np.max(self.meshy)
        print("    - X range [{}, {}]".format(x_min, x_max))
        print("    - Y range [{}, {}]".format(y_min, y_max))

        if self.ndim == 3:
            z_min = np.min(self.meshz)
            z_max = np.max(self.meshz)
            print("    - Z range [{}, {}]".format(z_min, z_max))

        if self.is_a_spectrum_file():
            print()
            print("  +> Spectrum information")
            print()

            freqs, _ = self.get_spectrum_freq()
            print("    - Number of frequencies: {}".format(len(freqs)))
            print("    - Frequencies range [{}, {}]"
                  .format(freqs.min(), freqs.max()))
            ndir = self.npoin3 - self.nelem3
            print("    - Number of directions : {} / angle:{} "
                  .format(ndir, 360./ndir))

        print()
        print("~> Parallel info")
        print()
        if self.nptir == 0:
            print("  - No parallel information")
        else:
            print("  - Number of interface points: {}".format(self.nptfr))
            print("  - Local ot Global numbering {}".format(self.knolg))

        print()
        print("~> Boundary info")
        print()
        if self.boundary_file is None:
            print("  - No boundary file given")
        else:
            print("  - Boundary element type: {}"
                  .format(elem2str(self.typ_bnd_elem)))

            print("  - Number of bnd points: {}".format(self.nptfr))
            print("  - Number of bnd elements: {}".format(self.nelebd))

            #TODO: display stuff on boundary values ?

        print()
        print("~> Data info")
        print()

        print("  - Number of records: {}".format(self.ntimestep))
        print("  - Time range: [{}, {}]".format(self.times[0], self.times[-1]))
        print("  - Number of variables: {}".format(self.nvar))

        for name, unit in zip(self.varnames, self.varunits):
            print("    - Name: {:16s} Unit: {}".format(name, unit))

        if full:
            print()
            for var in self.varnames:
                print("  - {}".format(var))
                for itime, time in enumerate(self.times):
                    val = self.get_data_value(var, itime)
                    print("    +> Time: {}s min, max:  [{}, {}]"
                          .format(time, val.min(), val.max()))
