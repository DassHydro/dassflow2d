# -*- coding: utf-8 -*-
# TODO: Add more logger info
"""
    Python wrapper to the Fortran APIs of module hermes of Telemac-Mascaret

    Author(s): Fabrice Zaoui, Yoann Audouin, Cedric Goeury, Renaud Barate

    Copyright EDF 2016
"""
import sys
import logging
from os import path
import numpy as np
from utils.exceptions import TelemacException
#
PRISM = 40
TETRAHEDRON = 30
QUADRANGLE = 20
TRIANGLE = 10
BND_SEGMENT = 55
BND_POINT = 1

def ndp2elem(ndp, ndim):
    """
    Identify element type from number of point per element

    @param ndp (int) Number of point per element
    @param ndim (int) Dimension of the element

    @returns ()
    """
    if ndp == 1:
        elem = BND_POINT
    elif ndp == 2:
        elem = BND_SEGMENT
    elif ndp == 3:
        elem = TRIANGLE
    elif ndp == 4:
        if ndim == 2:
            elem = QUADRANGLE
        else:
            elem = TETRAHEDRON
    elif ndp == 6:
        elem = PRISM
    else:
        elem = -1
    return elem


def elem2str(elem):
    """
    Return string version of elem from variable just before

    @param elem (int) Type of element

    @returns (str) Its name
    """
    string = ''
    if elem == PRISM:
        string = 'prism'
    elif elem == TETRAHEDRON:
        string = 'tetrahedron'
    elif elem == QUADRANGLE:
        string = 'quadrangle'
    elif elem == TRIANGLE:
        string = 'triangle'
    elif elem == BND_SEGMENT:
        string = 'bnd segment'
    elif elem == BND_POINT:
        string = 'bnd point'
    else:
        string = 'unknown'

    return string


class HermesFile():
    """The Generic Python class for TELEMAC-MASCARET APIs"""
    _hermes = None
    logger = logging.getLogger(__name__)
    _error = 0

    def __init__(self, file_name, fformat,
                 access='r', boundary_file=None,
                 log_lvl='INFO'):
        """
        Constructor for HermesFile

        @param file_name (string) Name of the file
        @param fformat (string) File format
        @param access (string) Access to the file ('r' for read 'w' for write)
        @param boundary_file (string) Name of the boundary file
        @param log_lvl (string) Logger level
        """

        if log_lvl == 'INFO':
            i_log = logging.INFO
        elif log_lvl == 'DEBUG':
            print("in debug")
            i_log = logging.DEBUG
        else:
            i_log = logging.CRITICAL
        logging.basicConfig(level=i_log)
        self.logger.setLevel(i_log)

        self.fformat = fformat.encode('utf-8') + b' '*(8 - len(fformat))
        self.file_name = file_name
        self.boundary_file = boundary_file

        try:
            self.logger.debug("Loading hermes f2py module")
            import _hermes
        except Exception as execpt:
            if sys.platform.startswith('linux'):
                ext = 'so'
            elif sys.platform.startswith('win'):
                ext = 'dll'
            else:
                raise TelemacException('Error: unsupported Operating System!')
            raise TelemacException(
                'Error: unable to load the dynamic library '
                + '_hermes.' + ext
                + '\nYou can check the environment variable:'
                + ' PYTHONPATH'
                + '\n'+str(execpt))
        HermesFile._hermes = sys.modules['_hermes']

        if 'r' in access:
            if 'w' in access:
                self.openmode = b'READWRITE'
            else:
                self.openmode = b'READ     '
                if not path.exists(self.file_name):
                    raise TelemacException(
                        "Could not find {}".format(self.file_name))
        elif 'w' in access:
            self.openmode = b'WRITE    '
            if path.exists(self.file_name):
                raise TelemacException(
                    "File already exist remove it first:{}"
                    .format(self.file_name))
        else:
            raise TelemacException(
                "Error in access string '%s' \
                should contain only r and/or w " % access)

        self.logger.debug("Opening mesh %s in format %s in mode %s",
                          self.file_name,
                          self.fformat,
                          self.openmode)
        self.my_id, self.error = \
            HermesFile._hermes.open_mesh(self.fformat,
                                         self.file_name,
                                         self.openmode)
        if self.boundary_file is not None:
            self.logger.debug("Opening bnd %s in format %s in mode %s",
                              self.fformat,
                              self.file_name,
                              self.openmode)
            self._errror = \
                HermesFile._hermes.open_bnd(self.fformat,
                                            self.boundary_file,
                                            self.my_id,
                                            self.openmode)

        if b'READ' in self.openmode:
            # Identifying elements type in files
            ndim = self.get_mesh_dimension()
            if ndim == 2:
                self.typ_elem = TRIANGLE
                # If no triangle it must be quadrangle
                nelem = self.get_mesh_nelem()
                if nelem == 0:
                    self.typ_elem = QUADRANGLE
            else:
                self.typ_elem = PRISM
            if 'med' in fformat.lower():
                self.typ_bnd_elem = BND_SEGMENT
            else:
                self.typ_bnd_elem = BND_POINT
        else:
            self.typ_elem = None
            self.typ_bnd_elem = None

        self.logger.debug("typ_elem: %s and typ_bnd_elem: %s",
                          elem2str(self.typ_elem),
                          elem2str(self.typ_bnd_elem))

    @property
    def error(self):
        """Error property
        """
        return self._error

    @error.setter
    def error(self, value):
        """Detect errors

        Overwright attribute setter to detect API errors.
        If :attr:`error` is not set null, an error is raised and the programme
        is terminated.

        :param int value: value to assign
        """
        if value != 0:
            self.logger.error("Hermes API error:\n%s",
                              HermesFile._hermes.get_error_message())
            raise TelemacException("Hermes API error:\n{}".format(
                HermesFile._hermes.get_error_message()))
        self._error = 0

    def close(self):
        """
        Closing file
        """
        if self.boundary_file is not None:
            self.logger.debug("Closing bnd file %s", self.boundary_file)
            if HermesFile._hermes is not None:
                self.error = HermesFile._hermes.close_bnd(self.fformat,
                                                          self.my_id)

        self.logger.debug("Closing mesh file %s", self.file_name)
        if HermesFile._hermes is not None:
            self.error = HermesFile._hermes.close_mesh(self.fformat,
                                                       self.my_id)

    def get_mesh_title(self):
        """
        Retuns the title of the file

        @returns The title
        """

        self.logger.debug("Getting title")
        title, self.error = HermesFile._hermes.get_mesh_title(self.fformat,
                                                              self.my_id)

        return title.decode('utf-8')

    def get_mesh_date(self):
        """
        Retuns the date of the file

        @returns The date (6-integer-array)
        """

        self.logger.debug("Getting date")
        date = np.zeros((6), dtype=np.int32)
        self.error = HermesFile._hermes.get_mesh_date(self.fformat,
                                                      self.my_id, date)

        return date

    def get_mesh_nelem(self):
        """
        Retuns the number of element in the file

        @returns The number of elements
        """

        self.logger.debug("Getting number of elements")
        nelem, self.error = \
            HermesFile._hermes.get_mesh_nelem(self.fformat,
                                              self.my_id,
                                              self.typ_elem)

        return nelem

    def get_mesh_npoin_per_element(self):
        """
        Retuns the number of points per element in the file

        @returns The number of points per element
        """

        self.logger.debug("Getting number of points per element")
        ndp, self.error = HermesFile._hermes.get_mesh_npoin_per_element(
            self.fformat, self.my_id, self.typ_elem)

        return ndp

    def get_mesh_connectivity(self):
        """
        Retuns the connectivity for the given type

        @returns An 2d array of shape (nelem, ndp)
        """

        self.logger.debug("Getting connectivity")
        nelem = self.get_mesh_nelem()
        self.logger.debug("Number of element: %d", nelem)
        ndp = self.get_mesh_npoin_per_element()
        self.logger.debug("Number of points per element: %d", ndp)
        tmp_ikle = np.zeros((nelem*ndp), dtype=np.int32)
        self.error = HermesFile._hermes.get_mesh_connectivity(
            self.fformat, self.my_id, self.typ_elem,
            tmp_ikle, nelem, ndp)
        ikle = tmp_ikle.reshape((nelem, ndp)) - 1

        return ikle

    def get_mesh_npoin(self):
        """
        Retuns the number of points

        @returns The number of points
        """

        self.logger.debug("Getting number of points %d %d",
                          self.my_id, self.typ_elem)
        npoin, self.error = HermesFile._hermes.get_mesh_npoin(
            self.fformat, self.my_id, self.typ_elem)

        return npoin

    def get_mesh_nplan(self):
        """
        Retuns the number of planes

        @returns The number of planes
        """

        self.logger.debug("Getting number of planes")
        nplan, self.error = HermesFile._hermes.get_mesh_nplan(
            self.fformat, self.my_id)

        return nplan

    def get_mesh_dimension(self):
        """
        Retuns the number of dimensions

        @returns The number of dimensions
        """

        self.logger.debug("Getting number of dimension")
        ndim, self.error = HermesFile._hermes.get_mesh_dimension(
            self.fformat, self.my_id)

        return ndim

    def get_mesh_orig(self):
        """
        Retuns the number of planes

        @returns The number of planes
        """

        self.logger.debug("Getting origin of coordinates")
        x_orig, y_orig, self.error = HermesFile._hermes.get_mesh_orig(
            self.fformat, self.my_id)

        return x_orig, y_orig

    def get_mesh_coord(self, jdim):
        """
        Retuns the coordinates of each points for a given dimension

        @param jdim Index of dimension [1-ndim]

        @returns A numpy array of size npoin
        """

        self.logger.debug("Getting coordinates for dimension %d", jdim)
        ndim = 0
        npoin = self.get_mesh_npoin()
        ndim = self.get_mesh_dimension()
        coord = np.zeros((npoin))
        # If in serafin and dimension 3 z coordinates is the
        # first variable at the first time step
        if jdim == 3 and b'SERAFIN' in self.fformat:
            var_names, _ = self.get_data_var_list()
            coord = self.get_data_value(var_names[0], 0)
        else:
            self.error = HermesFile._hermes.get_mesh_coord(
                self.fformat,
                self.my_id,
                jdim,
                ndim,
                coord,
                npoin)

        return coord

    def get_mesh_l2g_numbering(self):
        """
        Retuns the local to global numbering

        @returns The local to global numbering
        """

        self.logger.debug("Getting local to gloval numbering")
        npoin = self.get_mesh_npoin()
        knolg = np.zeros((npoin), dtype=np.int32)
        self.error = HermesFile._hermes.get_mesh_l2g_numbering(
            self.fformat,
            self.my_id,
            knolg,
            npoin)

        return knolg

    def get_mesh_nptir(self):
        """
        Retuns the number of interface points

        @returns The number of interface points
        """

        self.logger.debug("Getting number of interface points")
        nptir, self.error = HermesFile._hermes.get_mesh_nptir(
            self.fformat, self.my_id)

        return nptir

    def get_bnd_ipobo(self):
        """
        Retuns the ipobo array

        @returns The ipobo array
        """

        self.logger.debug("Getting boundary ipobo")
        npoin = self.get_mesh_npoin()
        nelebd = self.get_bnd_nelem()
        ipobo = np.zeros((npoin), dtype=np.int32)
        self.error = HermesFile._hermes.get_bnd_ipobo(
            self.fformat,
            self.my_id,
            nelebd,
            self.typ_bnd_elem,
            ipobo,
            npoin)

        return ipobo

    def get_bnd_numbering(self):
        """
        Retuns the boundary to general numbering

        @returns The boundary to general numbering
        """

        self.logger.debug("Getting boundary numbering")
        nptfr = self.get_bnd_npoin()
        nbor = np.zeros((nptfr), dtype=np.int32)
        self.error = HermesFile._hermes.get_bnd_numbering(
            self.fformat,
            self.my_id,
            self.typ_bnd_elem,
            nbor,
            nptfr)
        # Switching to Python Numbering
        nbor -= 1

        return nbor

    def get_bnd_connectivity(self):
        """
        Retuns the connectivity array for the boundary elements

        @returns The connectivity array for the boundary elements
        """
        self.logger.debug("Getting boundary connectivity")
        nelebd = self.get_bnd_nelem()
        if self.typ_bnd_elem == BND_SEGMENT:
            ndp = 2
        else:
            ndp = 1
        tmp_ikle_bnd = np.zeros((nelebd*ndp), dtype=np.int32)
        self.error = HermesFile._hermes.get_bnd_connectivity(
            self.fformat,
            self.my_id,
            self.typ_bnd_elem,
            nelebd, ndp,
            tmp_ikle_bnd)
        ikle_bnd = np.zeros((nelebd, ndp), dtype=np.int32)
        for i in range(nelebd):
            for j in range(ndp):
                ikle_bnd[i, j] = tmp_ikle_bnd[j*nelebd+i] -1
        del tmp_ikle_bnd

        return ikle_bnd

    def get_bnd_npoin(self):
        """
        Retuns the number of boundary points

        @returns The number of boundary points
        """

        self.logger.debug("Getting number of boundary points")
        nptfr, self.error = HermesFile._hermes.get_bnd_npoin(
            self.fformat,
            self.my_id,
            self.typ_bnd_elem)

        return nptfr

    def get_bnd_nelem(self):
        """
        Retuns the number of boundary elements

        @returns The number of boundary elements
        """

        self.logger.debug("Getting number of boundary elements")
        nelebd, self.error = HermesFile._hermes.get_bnd_nelem(
            self.fformat,
            self.my_id,
            self.typ_bnd_elem)

        return nelebd

    def get_bnd_value(self):
        """
        Retuns the information on the boundary values

        @returns liubor, lihbor, livbor, hbor, ubor, vbor, chbord,
                 litbor, tbor, atbor, btbor
        """

        self.logger.debug("Getting boundary values")
        nptfr = self.get_bnd_npoin()
        nelebd = self.get_bnd_nelem()
        lihbor = np.zeros((nptfr), dtype=np.int32)
        liubor = np.zeros((nptfr), dtype=np.int32)
        livbor = np.zeros((nptfr), dtype=np.int32)
        hbor = np.zeros((nptfr))
        ubor = np.zeros((nptfr))
        vbor = np.zeros((nptfr))
        chbord = np.zeros((nptfr))
        trac = True
        litbor = np.zeros((nptfr), dtype=np.int32)
        tbor = np.zeros((nptfr))
        atbor = np.zeros((nptfr))
        btbor = np.zeros((nptfr))
        self.error = HermesFile._hermes.get_bnd_value(
            self.fformat,
            self.my_id,
            self.typ_bnd_elem,
            nelebd,
            lihbor,
            liubor,
            livbor,
            hbor,
            ubor,
            vbor,
            chbord,
            trac,
            litbor,
            tbor,
            atbor,
            btbor,
            nptfr)

        color = np.zeros((nptfr), dtype=np.int32)
        self.error = HermesFile._hermes.get_bnd_color(
            self.fformat, self.my_id,
            self.typ_bnd_elem,
            color)

        return lihbor, liubor, livbor, hbor, ubor, vbor, chbord, \
            litbor, tbor, atbor, btbor, color

    def get_data_nvar(self):
        """
        Returns the number of variables

        @returns The number of variables
        """

        self.logger.debug("Get number of variables")
        nvar, self.error = HermesFile._hermes.get_data_nvar(
            self.fformat,
            self.my_id)

        return nvar

    def get_data_var_list(self):
        """
        Retuns the list of the variables name and units

        @returns Two arrays of size nvar first for name second for units
        """

        self.logger.debug("Get list of variable")
        nvar = self.get_data_nvar()
        res = HermesFile._hermes.get_data_var_list2(self.fformat,
                                                    self.my_id,
                                                    nvar)
        tmp_var_name, tmp_var_unit, self.error = res

        vnames = []
        vunit = []
        # Reordering string array for variable names
        # Extracting name and info into a list
        for i in range(nvar):
            var_name = b''.join(tmp_var_name[i*16:(i+1)*16]).decode('utf-8')
            var_unit = b''.join(tmp_var_unit[i*16:(i+1)*16]).decode('utf-8')
            vnames.append(var_name.strip())
            vunit.append(var_unit.strip())
        return vnames, vunit

    def get_data_ntimestep(self):
        """
        Retuns the number of time steps

        @returns The number of time steps
        """

        self.logger.debug("Get data ntimestep")
        ntimestep, self.error = HermesFile._hermes.get_data_ntimestep(
            self.fformat, self.my_id)

        return ntimestep

    def get_data_time(self, record):
        """
        Retuns the time of a given record

        @param record Number of the record (starts from 0)

        @returns The time
        """

        if record == -1:
            ntimestep = self.get_data_ntimestep()
            rrecord = ntimestep - 1
        else:
            rrecord = record

        self.logger.debug("Get data time at %d", record)
        time, self.error = HermesFile._hermes.get_data_time(
            self.fformat,
            self.my_id,
            rrecord)

        return time

    def get_data_value(self, var_name, record):
        """
        Retuns the value for each point for a given variable and a given record

        @param var_name Name of the variable
        @param record Number of the record (starts from 0 if a -1 is given will
                                   give the last time step)

        @returns A numpy array of size npoin
        """

        npoin = self.get_mesh_npoin()
        values = np.zeros((npoin))
        var_name2 = var_name + ' '*(16 - len(var_name))

        if record == -1:
            ntimestep = self.get_data_ntimestep()
            rrecord = ntimestep - 1
        else:
            rrecord = record

        self.logger.debug("Getting data for %s at %d",
                          var_name2, record)
        self.error = HermesFile._hermes.get_data_value(
            self.fformat,
            self.my_id,
            rrecord,
            var_name2,
            values,
            npoin)

        return values

    def set_header(self, title, nvar, var_name, var_unit):
        """
        Write header of the file

        @param title Title of the file
        @param nvar Number of variables
        @param var_name Name for each variable
        @param var_unit Unit for each variable
        """
        if len(title) > 80:
            tmp_title = title[:80].encode('utf-8')
        else:
            tmp_title = title.encode('utf-8') + b' '*(80-len(title))
        tmp_var_name = [b' ']*32*nvar
        for i, (var, unit) in enumerate(zip(var_name, var_unit)):
            for j, varj in enumerate(var):
                tmp_var_name[i*32+j] = varj.encode('utf-8')
            for j, unitj in enumerate(unit):
                tmp_var_name[i*32+16+j] = unitj.encode('utf-8')

        self.logger.debug("Writing header information")
        self.error = HermesFile._hermes.set_header(
            self.fformat, self.my_id,
            tmp_title, tmp_var_name, nvar)

    def set_mesh(self, mesh_dim, typ_elem, ndp, nptfr, nptir, nelem, npoin,
                 ikles, ipobo, knolg, coordx, coordy, nplan, date,
                 time, x_orig, y_orig, coordz=None):
        """
        Write the mesh information into the file

        @param mesh_dim
        @param mesh_dim Dimension of the mesh
        @param typ_elem TYPE OF THE MESH ELEMENTS
        @param ndp Number of points per element
        @param nptfr Number of boundary point
        @param nptir Number of interface point
        @param nelem Number of element in the mesh
        @param npoin Number of points in the mesh
        @param ikles Connectivity array for the main element
        @param ipobo Is a boundary point ? array
        @param knolg Local to global numbering array
        @param coordx X coordinates of the mesh points
        @param coordy Y coordinates of the mesh points
        @param nplan Number of planes
        @param date Date of the creation of the mesh
        @param time Time of the creation of the mesh
        @param x_orig Origin of x coordinates
        @param y_orig Origin of y coordinates
        @param coordz Z coordinates of the mesh points
        """
        if coordz is None:
            tmp_z = np.zeros((npoin), dtype=np.float64)
        else:
            tmp_z = coordz

        self.typ_elem = typ_elem

        tmp_ikle = ikles.T.reshape((nelem*ndp)) + 1

        self.logger.debug("Writing mesh information")
        self.error = HermesFile._hermes.set_mesh(
            self.fformat, self.my_id,
            mesh_dim, typ_elem, ndp, nptfr,
            nptir, nelem, tmp_ikle,
            ipobo, knolg, coordx, coordy,
            nplan, date, time, x_orig, y_orig,
            npoin, tmp_z)
        del tmp_ikle
        if coordz is None:
            del tmp_z

    def add_data(self, var_name, var_unit, time, record, first_var, values):
        """
        Write inform ation for a given variable and a given timestep

        @param var_name Name of the variable
        @param var_unit Unit of the variable
        @param time Time of the data
        @param record Time step of the data (starts from 0)
        @param first_var True if it is the first variable of the dataset
        @param values The value for each point of the mesh
        """
        nval = len(values)
        tmp_var_name = var_name.encode('utf-8') + b' '*(16-len(var_name)) +\
            var_unit.encode('utf-8') + b' '*(16-len(var_unit))

        self.logger.debug("Writing data for %s at record %d",
                          tmp_var_name, record)
        self.error = HermesFile._hermes.add_data(
            self.fformat, self.my_id,
            tmp_var_name, time, record,
            first_var, values, nval)

    def set_bnd(self, typ_bnd_elem, nelebd, ikle, lihbor, liubor,
                livbor, hbor, ubor, vbor, chbord, litbor, tbor, atbor,
                btbor, color):
        """
        Write boundary information

        @param typ_bnd_elem (int) Type of boundary element
        @param nelebd Number of boundary elements
        @param ikle Connectivity array for the boundary elements
        @param lihbor Type of boundary conditions on depth
        @param liubor Type of boundary conditions on u
        @param livbor Type of boundary conditions on v
        @param hbor Prescribed boundary condition on depth
        @param ubor Prescribed boundary condition on velocity u
        @param vbor Prescribed boundary condition on velocity v
        @param chbord Friction coefficient at boundary
        @param litbor Physical boundary conditions for tracers
        @param tbor Prescribed boundary condition on tracer
        @param atbor Thermal exchange coefficients
        @param btbor Thermal exchange coefficients
        @param color Boundary color of the boundary element
        """
        self.typ_bnd_elem = typ_bnd_elem
        if self.typ_bnd_elem == BND_SEGMENT:
            ndp = 2
        else:
            ndp = 1
        # Switching
        tmp_ikle = ikle.T.reshape((nelebd*ndp)) + 1

        self.logger.debug("Writing boundary file")
        self.error = HermesFile._hermes.set_bnd(
            self.fformat, self.my_id,
            typ_bnd_elem, nelebd, ndp, tmp_ikle,
            lihbor,
            liubor, livbor, hbor, ubor, vbor,
            chbord, litbor, tbor, atbor, btbor,
            color)

        del tmp_ikle


    def import_group_info(self, src):
        """
        Import group information from a given HermesFile

        @param src HermesFile from which group will be imported
        """
        self.typ_bnd_elem = src.typ_bnd_elem
        if self.typ_bnd_elem is None:
            tmp_typ_bnd_elem = 0
        else:
            tmp_typ_bnd_elem = self.typ_bnd_elem
        ikle_bnd = src.get_bnd_connectivity()
        nelebd, ndp = ikle_bnd.shape
        tmp_ikle_bnd = ikle_bnd.T.reshape((nelebd*ndp))

        self.logger.debug("Transfering group information")
        self.error = HermesFile._hermes.transfer_group_info(
            self.fformat, src.my_id,
            self.my_id, self.typ_elem,
            tmp_typ_bnd_elem, tmp_ikle_bnd, nelebd, ndp,
            False, False)


    def import_from(self, src):
        """
        Rewriting src into dst

        @param src HermesFile from which the data is imported
        """
        self.logger.debug("Getting information from file %s", src.file_name)
        title = src.get_mesh_title()
        date = src.get_mesh_date()

        ndim = src.get_mesh_dimension()
        typ_elem = src.typ_elem
        typ_bnd_elem = src.typ_bnd_elem

        npoin = src.get_mesh_npoin()
        nelem = src.get_mesh_nelem()
        ndp = src.get_mesh_npoin_per_element()
        nplan = src.get_mesh_nplan()
        x_orig, y_orig = src.get_mesh_orig()

        coordx = src.get_mesh_coord(1)
        coordy = src.get_mesh_coord(2)
        coordz = None
        if ndim == 3:
            coordz = src.get_mesh_coord(3)

        ikle = src.get_mesh_connectivity()

        nptir = src.get_mesh_nptir()
        if nptir > 0:
            knolg = src.get_mesh_l2g_numbering()
        else:
            knolg = np.zeros((npoin), dtype=np.int32)

        if src.boundary_file is not None:
            typ_bnd_elem = src.typ_bnd_elem
            nptfr = src.get_bnd_npoin()
            nelebd = src.get_bnd_nelem()
            ipobo = src.get_bnd_ipobo()
            ikle_bnd = src.get_bnd_connectivity()

            lihbor, liubor, livbor, hbor, ubor, vbor, chbord, \
                litbor, tbor, atbor, btbor, color = src.get_bnd_value()
        else:
            nptfr = 0
            nelebd = 0
            ipobo = np.zeros((npoin), dtype=np.int32)

        ntimestep = src.get_data_ntimestep()
        nvar = src.get_data_nvar()
        var_name, var_unit = src.get_data_var_list()

        self.set_header(title, nvar, var_name, var_unit)

        date2 = np.zeros((3), dtype=np.int32)
        time2 = np.zeros((3), dtype=np.int32)
        date2[0] = date[0]
        date2[1] = date[1]
        date2[2] = date[2]
        time2[0] = date[3]
        time2[1] = date[4]
        time2[2] = date[5]

        self.set_mesh(ndim, typ_elem, ndp, nptfr, nptir, nelem, npoin,
                      ikle, ipobo, knolg, coordx, coordy, nplan, date2,
                      time2, x_orig, y_orig, coordz)

        if self.fformat == src.fformat and b'SERAFIN' not in self.fformat:
            self.import_group_info(src)
        else:
            if src.boundary_file is not None:
                self.set_bnd(typ_bnd_elem, nelebd, ikle_bnd, lihbor, liubor,
                             livbor, hbor, ubor, vbor, chbord, litbor, tbor,
                             atbor, btbor, color)

        for i in range(ntimestep):
            time = src.get_data_time(i)
            for j in range(nvar):
                values = src.get_data_value(var_name[j], i)
                self.add_data(var_name[j], var_unit[j],
                              time, i, j == 0, values)
                del values

    def __repr__(self):
        """
        representation of the object
        """
        string = '*'*32 + '\n'
        string += 'Generic info' + '\n'
        string += '*'*32 + '\n'

        string += "Title: %s\n" % self.get_mesh_title()[1:72]
        date = self.get_mesh_date()
        string += "Date: %d/%d/%d %dH%dM%dS\n" % (date[2], date[1], date[0],
                                                  date[3], date[4], date[5])

        string += '*'*32 + '\n'
        string += 'Mesh info\n'
        string += '*'*32 + '\n'

        ndim = self.get_mesh_dimension()
        string += "Ndim: %d\n" % ndim

        if self.typ_elem == TRIANGLE:
            string += "Element type: TRIANGLE\n"
        elif self.typ_elem == PRISM:
            string += "Element type: PRISM\n"
        else:
            string += "Element type: UNKNOWN\n"
        string += "Npoin: %d\n" % self.get_mesh_npoin()
        string += "Nelem: %d\n" % self.get_mesh_nelem()
        string += "Ndp: %d\n" % self.get_mesh_npoin_per_element()
        string += "nplan: %d\n" % self.get_mesh_nplan()

        string += "coordinates:\n"
        string += " - On x :%s\n" % str(self.get_mesh_coord(1))
        string += " - On y :%s\n" % str(self.get_mesh_coord(2))
        if ndim == 3:
            string += " - On z :%s\n" % str(self.get_mesh_coord(3))

        string += "ikle: %s\n" % str(self.get_mesh_connectivity())

        string += '*'*32 + '\n'
        string += 'Parallel info\n'
        string += '*'*32 + '\n'
        nptir = self.get_mesh_nptir()
        if nptir != 0:
            string += "Nptir: %d\n" % nptir
            string += "knolg: %s\n" % str(self.get_mesh_l2g_numbering())
        else:
            string += 'No parallel information\n'

        string += '*'*32 + '\n'
        string += 'Bnd info\n'
        string += '*'*32 + '\n'
        if self.boundary_file is not None:
            if self.typ_bnd_elem == BND_POINT:
                string += "Bnd element: BND_POINT\n"
            elif self.typ_bnd_elem == BND_SEGMENT:
                string += "Bnd element: BND_SEGMENT\n"
            else:
                string += "Bnd element: UNKNOWN\n"
            string += "Nptfr: %d\n" % self.get_bnd_npoin()
            string += "Nelebd: %d\n" % self.get_bnd_nelem()

            string += "ipobo: %s\n" % self.get_bnd_ipobo()
            string += "nbor: %s\n" % self.get_bnd_numbering()
            string += "ikle_bnd: %s\n" % self.get_bnd_connectivity()

            lihbor, liubor, livbor, _, _, _, _, \
                litbor, _, _, _, _ = self.get_bnd_value()
            string += "bnd_values: \n"
            string += str(zip(lihbor, liubor, livbor, litbor)) + '\n'
        else:
            string += "No boundary information"

        string += '*'*32 + '\n'
        string += 'Data info\n'
        string += '*'*32 + '\n'

        ntimestep = self.get_data_ntimestep()
        string += "ntimestep: %d\n" % ntimestep
        nvar = self.get_data_nvar()
        string += "nvar: %d\n" % nvar
        var_name, var_unit = self.get_data_var_list()
        string += "var info:\n"
        for var in zip(var_name, var_unit):
            string += str(var) + '\n'

        for i in range(ntimestep):
            string += "Time: %fs\n" % self.get_data_time(i)
            for j in range(nvar):
                string += " - for %s:\n" % var_name[j]
                string += "     " + str(self.get_data_value(var_name[j], i))\
                          + '\n'

        return string
