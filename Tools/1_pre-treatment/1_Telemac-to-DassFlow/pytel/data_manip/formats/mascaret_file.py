"""
Contains classe to parse mascaret files
"""
from collections import OrderedDict
import csv
import logging
import struct
import os
from math import sqrt
import traceback
import sys

from utils.exceptions import MascaretException
import numpy as np


class Reach:
    """
    Reach object containing sections

    id (int) reach identifier
    name (str) reach name
    section ([Section]) list of sections
    """

    def __init__(self, my_id, name=None):
        self.id = my_id
        self.name = name
        if self.name is None:
            self.name = 'Reach_%i' % my_id
        self.sections = OrderedDict()
        self.nsections = 0
        self._n = None

    def add_section(self, section):
        """
        Add section
        @param section (Section)
        @return:
        """
        self.nsections += 1
        self.sections[section.id] = section

    def get_section_pk_list(self):
        """
        Get list of pk section
        @return (numpy.array)
        """
        return np.array([section.pk for _, section in self.sections.items()])

    def get_section_id_list(self):
        """
        Get list of id section
        @return (numpy.array)
        """
        return np.array([section.id for _, section in self.sections.items()])

    def get_section_idx(self, section_id):
        """
        Get section index
        @param section_id (int) section identifier
        @return (int) request section index
        """
        section_ids = list(self.sections.keys())
        try:
            return section_ids.index(section_id)
        except ValueError:
            raise MascaretException('Section identifier %i is not found.\n'
                                    'Possible section identifiers are:\n%s'
                                    % (section_id, section_ids))

    def __repr__(self):
        return 'Reach #%i (%s) with %i sections' % \
                (self.id, self.name, self.nsections)

    # Enable iteration over section list
    def __iter__(self):
        self._n = 0
        return self

    def __next__(self):
        if self._n > self.nsections - 1:
            raise StopIteration
        result = self.__getitem__(self._n)
        self._n += 1
        return result

    def __getitem__(self, i):
        if isinstance(i, slice):
            start = 0 if i.start is None else i.start
            stop = -1 if i.stop is None else i.stop
            step = 1 if i.step is None else i.step
            return [self[j] for j in range(start, stop, step)]
        return list(self.sections.items())[i][1]


class Section:
    """
    Geometry of a cross-section

    id (int) profile identifier
    name (str) profile name
    pk (float) distance along the hydraulic axis

    axis (numpy 1D-array) coordinates of hydraulic axis
    x (numpy 1D-array) point coordinates along x axis
    y (numpy 1D-array) point coordinates along y axis
    distances (numpy 1D-array) cumulative distance from first point along the
    profile
    nb_points (int) number of points
    limits ({limit_name: point_numbering}) position of limits
    """

    def __init__(self, my_id, pk, name=None):
        self.id = my_id
        self.name = name
        if self.name is None:
            self.name = 'Profil_%i' % my_id
        self.pk = pk

        self.axis = None
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.distances = np.array([])
        self.nb_points = 0
        self.limits = {}

        self.nlayers = 0
        self.layers_elev = None
        self.layer_names = []

    def set_points_from_trans(self, dist_array, z_array):
        if len(dist_array) != len(z_array):
            raise MascaretException(
                 'set_points_from_trans:\
                  Input arrays have not the same length')
        self.allocate(len(dist_array))
        for i, (dist, z) in enumerate(zip(dist_array, z_array)):
            if i == 0:
                limit = 'RD'
            elif i == (self.nb_points - 1):
                limit = 'RG'
            else:
                limit = None
            self.set_point(i, self.pk, dist, z, limit)

    def set_points_from_xyz(self, x_list, y_list, z_list):
        if not len(x_list) == len(y_list) == len(z_list):
            raise MascaretException(
                'set_points_from_xyz: Input arrays have not the same length')
        self.allocate(len(x_list))
        for i, (x, y, z) in enumerate(zip(x_list, y_list, z_list)):
            if i == 0:
                limit = 'RD'
            elif i == (self.nb_points - 1):
                limit = 'RG'
            else:
                limit = None
            self.set_point(i, x, y, z, limit)

    def set_axis(self, xa, ya):
        self.axis = (xa, ya)

    def get_limit(self, limit_name):
        try:
            return self.limits[limit_name]
        except KeyError:
            raise MascaretException('Limit %s is not found in %s' %
                                    (limit_name, self))

    def point_index_limit(self, i):
        for limit_name, index in self.limits.items():
            if index == i:
                return limit_name
        return None

    def allocate(self, nb_points):
        self.x = np.empty(nb_points)
        self.y = np.empty(nb_points)
        self.z = np.empty(nb_points)
        self.distances = np.empty(nb_points)
        self.nb_points = nb_points

    def set_point(self, i, x, y, z, limit=None):
        self.x[i] = x
        self.y[i] = y
        self.z[i] = z
        if limit is not None:
            self.limits[limit] = i
        if i == 0:
            self.distances[i] = 0
        else:
            self.distances[i] = self.distances[i - 1] + \
                                sqrt((self.x[i] - self.x[i - 1]) ** 2 +
                                     (self.y[i] - self.y[i - 1]) ** 2)

    def add_layer(self, thickness_table, name=None):
        self.nlayers += 1
        if name is None:
            name = 'Layer ' + str(self.nlayers)
        self.layer_names.append(name)
        if self.layers_elev is None:
            self.layers_elev = np.empty((1, self.nb_points))
            self.layers_elev[0, :] = self.z - thickness_table

        else:
            self.layers_elev = np.vstack((self.layers_elev,
                                          self.layers_elev[self.nlayers - 2]
                                          - thickness_table))

    def iter_on_points(self):
        for i, (x, y, z) in enumerate(zip(self.x, self.y, self.z)):
            limit = self.point_index_limit(i)
            limit_str = limit if limit is not None else ''
            yield x, y, z, limit_str

    def common_limits(self, other):
        """
        :param other: upstream or downstream section
        :type other: Section
        :return: list of limits
        """
        return list(set(self.limits.keys()).intersection(other.limits.keys()))

    def check_elevations(self):
        pass  # TODO

    def __repr__(self):
        return 'Section #%i (%s) at pk %f' % (self.id, self.name, self.pk)


class MascaretFileParent:
    """
    Parse MascaretFile ('opt' or 'rub')
    Warning: Encoding is not checked and depends on system configuration
    """

    logger = logging.getLogger(__name__)

    def __init__(self, file_name, access='r', log_lvl='INFO'):
        """
        Constructor for MascaretFile
        Warning Only suited for results at cross-sections (not adapted to
        Casier or Traceur outputs)

        @param file_name (str) Name of the file
        @param access (str) Access to the file ('r' for read 'w' for write, add
        'b' for binary file)
        @param log_lvl (str) logger level

        Attributs:
        - file_name: file name
        - nreach: number of reaches
        - nvar: number of variables
        - varnames: variable names
        - varunits: unit of variables
        - times: list of time in seconds
        - _times_pos: starting position of frame in file
        """
        if log_lvl == 'INFO':
            i_log = logging.INFO
        elif log_lvl == 'DEBUG':
            i_log = logging.DEBUG
        else:
            i_log = logging.CRITICAL
        logging.basicConfig(level=i_log)

        # File name
        self.file_name = file_name
        self._position_first_frame = 0

        # Attributes for geometry
        self.nreaches = 0
        self.nsections = 0
        self._reaches = OrderedDict()

        # Attributes for variables
        self.nvar = 0
        self.varnames_dict = {'names': [],
                              'abbr': [],
                              'id': [],
                              'units': []}

        # Attributes for section variables
        self.nsectionvar = 0
        self.section_varnames_dict = {'names': [],
                                      'abbr': [],
                                      'id': [],
                                      'units': []}

        # Attributes for temporal data
        self._times = []
        self._times_pos = []
        self._ntimestep = 0

        self._file = open(self.file_name, access)

    def __del__(self):
        self.logger.debug("Closing mesh file %s", self.file_name)
        if self._file is not None:
            self._file.close()

    def __repr__(self):
        return 'MascaretFile: %s (mode=%s)' % (self.file_name, self._file.mode)

    @property
    def times(self):
        """
        Returns a list of the times in the file
        """
        if not self._times:
            self.get_time()
        return self._times

    @property
    def times_pos(self):
        """
        Returns a list of the frame position in the file
        """
        if not self._times_pos:
            self.get_time()
        return self._times_pos

    @property
    def ntimestep(self):
        """
        Returns a list of the times step number in the file
        """
        if self._ntimestep == 0:
            self.get_time()
        return self._ntimestep

    @property
    def reaches(self):
        """
        Returns a list of the reaches in the file
        """
        if not self._reaches:
            self.get_reaches()
        return self._reaches

    def _move_to_first_frame(self):
        """Start file reader position before first frame
         definition for Rubens file"""
        self._file.seek(self._position_first_frame)

    def get_time(self):
        raise NotImplementedError('Has to be override in subclass')

    def get_reaches(self):
        raise NotImplementedError('Has to be override in subclass')

    def error(self, message):
        raise MascaretException('ERROR: %s' % message)

    def add_variable(self, varname, varunit, varname_abbr):
        """Add a variable"""
        self.varnames_dict['names'].append(varname)
        self.varnames_dict['abbr'].append(varname_abbr)
        self.varnames_dict['units'].append(varunit)
        self.varnames_dict['id'].append(self.nvar)
        self.nvar += 1

    def add_section_variable(self, varname, varunit, varname_abbr):
        """Add a section variable"""
        self.section_varnames_dict['names'].append(varname)
        self.section_varnames_dict['abbr'].append(varname_abbr)
        self.section_varnames_dict['units'].append(varunit)
        self.section_varnames_dict['id'].append(self.nsectionvar)
        self.nsectionvar += 1

    def get_position_var(self, var_name, type='names'):
        """
        Get position variable

        @param var_name (string) variable name
        @param type (string) 'names' for full name or 'abbr' for abreviation

        @return variable name index
        """
        try:
            if type == 'abbr':
                return self.varnames_dict['abbr'].index(var_name)
            else:
                return self.varnames_dict['names'].index(var_name)
        except ValueError:
            self.error('Variable `%s` not found. '
                       'Possibles values are:\n%s' %
                       (var_name, self.varnames_dict[type]))

    def get_values_at_reach(self, record, reach_id, vars_indexes=None):
        """
        Get values for variables for a single reach
        @param record (int) time index
        @param reach_id (int) reach index
        @param vars_indexes (list) List of variable names
        @return (numpy.array)
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        return self.get_values(record, vars_indexes)[reach_id]

    def get_values_at_sections(self, record, section_id, reach_id=1,
                               section_vars_indexes=None):
        """
        Get values for section variables for a single reach

        @param record (int) time index
        @param section_id (int) id of the section
        @param reach_id (int) reach index
        @param section_vars_indexes (list) List of variable names

        @return (numpy.array)
        """

        if self.fformat != 'ptravers':
            self.error("The method get_values_at_section is not possible\
                        with the kind of file loaded"
                       " (only for ptravers which stores section variables)")

        if self.nsectionvar == 0:
            self.error('There is no section variable loaded yet')

        if section_vars_indexes is None:
            section_vars_indexes = self.section_varnames_dict['id']
        _, res_sections = \
            self.get_values(
                            record,
                            get_section_values=True,
                            section_vars_indexes=section_vars_indexes
                           )
        return res_sections[reach_id][section_id]

    def get_position_var_abbr(self, var_abbr):
        """
        Get position  of the abbreviation variable

        @param var_abbr (string) variable name

        @return variable name index
        """
        return self.varnames_dict['abbr'].index(var_abbr)

    def get_reach_id(self, reach_name):
        """
        Get id position reach
        @param reach_name (string) reach name
        @return: reach name index
        """
        for reach in self.reaches:
            if reach_name == reach.name:
                return reach.id
        raise MascaretException('Reach name not found')

# TODO: add section variables (see write_opt_file)
    def write_optfile_header(self, outfile, vars_indexes=None):
        """
        write header file

        @param outfile (File) file open object
        @param vars_indexes (list) List of var index to write If none given
        writes them all
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']
        outfile.write('[variables]\n')
        for i in vars_indexes:
            outfile.write('"{0}";"{1}";"{2}";0\n'.format(
                    self.varnames_dict['names'][i],
                    self.varnames_dict['abbr'][i],
                    self.varnames_dict['units'][i])
                         )
        outfile.write('[resultats]\n')

# TODO: add section variables (mean value for each section?)
    def write_optfile_frame(self, outfile, res, time):
        """
        Write a one frame in opthyca file

        @param outfile (string) Name of output file
        @param res (OrderedDict) results data
        @param time (float) time value
        @return
        """
        for key in res.keys():
            id = self.reaches[key].get_section_id_list()
            pk = self.reaches[key].get_section_pk_list()
            for id, pk, val in zip(id, pk, res[key]):
                outfile.write(
                    '{0};"{1:2}";"{2:5}";{3};{4} \n'
                    .format(time, key, id, pk, ";"
                            .join([str(var) for var in val])))

# TODO: add an option to write also section variables
# (mean value for each section?)
    def write_optfile(self, outfile_name, times_indexes=None,
                      vars_indexes=None, timecheck=False):
        """
        Write an output file in opthyca format

        @param outfile_name (string) output file name
        @param times_indexes (list) List of time step index
        @param vars_indexes (list) List of variable names
        @param timecheck (boolean) Tf True more verbose
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']
        if times_indexes is None:
            times_indexes = [i for i in range(len(self.times))]

        outfile = open(outfile_name, 'w')
        # write header
        self.write_optfile_header(outfile, vars_indexes)
        # write times
        for id_time in times_indexes:
            res = self.get_values(id_time, vars_indexes)
            # write time step
            if timecheck:
                print('Writing timestep: ', self.times[id_time])
            self.write_optfile_frame(outfile, res, self.times[id_time])

        outfile.close()

    def export_as_lig(self, file_name, record):
        """
        Write Mascaret restart file

        @param file_name: (string) file name
        @param record (int) time index
        @return
        """
        pk = np.array([])
        i1 = {}
        i2 = {}
        nsection = 0
        for _, reach in self.reaches.items():
            pk = np.concatenate((pk, reach.get_section_pk_list()))
            i1[reach.id] = min([section.id for section in reach])
            i2[reach.id] = max([section.id for section in reach])
            nsection += reach.nsections

        i1i2 = []
        for b in sorted(i1.keys()):
            i1i2.append(str(i1[b]))
            i1i2.append(str(i2[b]))

        z_pos = self.get_position_var_abbr('Z')
        q_pos = self.get_position_var_abbr('Q')
        res = self.get_values(record, [z_pos, q_pos])
        z = np.array([])
        q = np.array([])
        for k in res.keys():
            z = np.concatenate((z, res[k][:, 0]))
            q = np.concatenate((q, res[k][:, 1]))

        result = {}
        result['X'] = pk
        result['Z'] = z
        result['Q'] = q

        with open(file_name, 'w') as fich:
            # Date is hardcoded, but it could be defined as:
            # datetime.datetime.utcnow()
            fich.write(
                'RESULTATS CALCUL,DATE : 01/01/1900 00:00\n')
            fich.write('FICHIER RESULTAT MASCARET{0}\n'.format(' ' * 47))
            fich.write('{0} \n'.format('-' * 71))
            fich.write(' IMAX  = {0:4} NBBIEF= {1:3}\n'
                       .format(str(nsection), str(self.nreaches)))
            chaine = [""]
            for k in range(0, len(i1i2), 10):
                chaine.append('I1,I2 =')
                for i in range(k, k + 10):
                    if i < len(i1i2):
                        chaine.append('{0:4}'.format(i1i2[i]))
                chaine.append("\n")
            fich.write(" ".join(chaine))

            for k in ['X', 'Z', 'Q']:
                fich.write(' ' + k + '\n')
                long = 0
                for x in result[k]:
                    fich.write('{:13.6f}'.format(x))
                    long += 1
                    if long == 5:
                        fich.write('\n')
                        long = 0

                if long != 0:
                    fich.write('\n')

            fich.write(' FIN\n')

# TODO: add section variables in the summary
    def summary(self):
        txt = '~> %s\n' % self
        for _, reach in self.reaches.items():
            txt += '    - %s\n' % reach
            for i, section in enumerate(reach):
                txt += '        %i) %s\n' % (i, section)
        txt += '%i variables:\n' % self.nvar
        for i, varname in enumerate(self.varnames_dict['names']):
            txt += '    - %s (%s) | %s \n' % (varname,
                                              self.varnames_dict['units'][i],
                                              self.varnames_dict['abbr'][i]
                                              )
        if self.nsectionvar != 0:
            txt += '%i section variables:\n' % self.nsectionvar
            for i, varname in enumerate(self.section_varnames_dict['names']):
                txt += '    - %s (%s) | %s \n' %\
                    (varname,
                     self.section_varnames_dict['units'][i],
                     self.section_varnames_dict['abbr'][i])
        txt += '%i temporal frames:\n' % self.ntimestep
        for i, time in enumerate(self.times):
            txt += '    - %i) %f\n' % (i, time)
        return txt


class Opthyca(MascaretFileParent):

    def __init__(self, file_name, access='r', log_lvl='INFO'):
        """
        Constructor for Opthyca file
        Warning: Only suited for results at cross-sections (not adapted to
        Casier or Traceur outputs)

        @param file_name (string) Name of the file
        @param access (string) Access to the file ('r' for read 'w' for write)
        @param log_lvl (string) Logger level

        Attributs specified to Opthyca:
        - fformat
        """
        super().__init__(file_name, access=access, log_lvl=log_lvl)
        self.fformat = 'opt'

        self._read_variables()

    def read_line(self):
        return self._file.readline().rstrip('\n')

    def _read_variables(self):
        """Read header with variables"""
        # Skip comments before variable definition
        row = self.read_line()
        while row != '[variables]':
            row = self.read_line()

        # Read variable definitions
        row = self.read_line()
        while row != '[resultats]':
            try:
                name, abbr, unit, _ = row.split(';')
            except ValueError:
                self.error('Variable description is not readable')
            self.add_variable(name.strip('\"'), unit.strip('\"'),
                              abbr.strip('\"'))
            row = self.read_line()
        self._position_first_frame = self._file.tell()

    def _read_line_resultat(self):
        """Interpret a line containing some results"""
        row = self.read_line()

        # To handle the case where there is only one time step
        if row == '':
            return np.nan, 'no_bief', np.nan, 'no_data'

        try:
            time_str, bief_name, _, pk_str, values_str = \
                    row.split(';', maxsplit=4)
        except ValueError:
            self.error('Number of values (separated by a semi-colon)\
             has to be more than 4!')

        try:
            time = float(time_str)
            section_pk = float(pk_str)
            values = [float(x) for x in values_str.split(';')]
        except ValueError as e:
            self.error(str(e))
        if len(values) != self.nvar:
            self.error('Number of values not coherent:\
             %i instead of %i' % (len(values), self.nvar))

        return time, int(bief_name.strip().strip('\"')), section_pk, values

    def get_reaches(self):
        """Read geometry"""
        self._move_to_first_frame()
        self.nreaches = 0

        time, reach_id, pk, _ = self._read_line_resultat()
        reach = Reach(reach_id)
        section_id = 1
        first_time = time
        prev_bief_id = reach_id

        while time == first_time:
            if prev_bief_id != reach_id:
                self.nreaches += 1
                self._reaches[reach.id] = reach
                reach = Reach(reach_id)

            prev_bief_id = reach_id
            reach.add_section(Section(section_id, pk))
            section_id += 1
            time, reach_id, pk, _ = self._read_line_resultat()
        self.nreaches += 1
        self._reaches[reach.id] = reach

    def get_values(self, record, vars_indexes=None):
        """
        Get values for variables for a given time index

        @param record (int) time index
        @param vars_indexes (list) List of variable names
        @return (dict) dict of 2D array with reach.id as key
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        res = OrderedDict()
        requested_time = self.times[record]
        self._file.seek(self._times_pos[record])
        for _, reach in self.reaches.items():
            all_values = []
            for _, section in reach.sections.items():
                time, reach_id, pk, values = self._read_line_resultat()
                if reach_id != reach.id:
                    self.error("Reach identifier is not coherent")
                if time != requested_time:
                    self.error("Time is not not coherent")
                if section.pk != pk:
                    self.error("pk is not not coherent")
                all_values.append(np.array(values)[vars_indexes])

            res[reach.id] = np.asarray(all_values)
        return res

    def get_series(self, reach_id, section_id, vars_indexes=None):
        """
        Get values for all variables for a give reach index and a given section
        index
        @param reach_id (int) reach index
        @param section_id (int) section index
        @param vars_indexes (list) List of variable names
        @return (numpy.array)
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        all_values = []
        section_idx = self.reaches[reach_id].get_section_idx(section_id)
        for record in range(self.ntimestep):
            tmp = self.get_values(record)
            all_values.append(tmp[reach_id][section_idx, vars_indexes])

        return np.array(all_values)

# TODO: verify this is working even if reaches has not been initialized
    def get_time(self):
        """
        Initialize time variables
        """
        self._move_to_first_frame()
        pos = self._file.tell()

        line = self.read_line()

        while line != '':
            time = float(line.split(';')[0])
            self._times.append(time)
            self._times_pos.append(pos)
            self._ntimestep += 1
            for _, reach in self.reaches.items():
                for _ in range(reach.nsections):
                    pos = self._file.tell()
                    line = self.read_line()


class Rubens(MascaretFileParent):

    def __init__(self, file_name, access='r', log_lvl='INFO'):
        """
        Constructor for Rubens file

        @param file_name (string) Name of the file
        @param access (string) Access to the file ('r' for read 'w' for write)
        @param log_lvl (string) Logger level

        Attributs specific to Rubens:
        - fformat
        - _size_file
        - _endians
        """
        super().__init__(file_name, access=access + 'b', log_lvl=log_lvl)
        self.fformat = 'rub'
        self._size_file = os.path.getsize(self.file_name)
        self._endians = ''

        self._read_binary_header()

    @property
    def endians(self):
        """
        Returns the endians of the file
        """
        if not self._endians:
            self.get_endians()
        return self._endians

    @staticmethod
    def _read_dico_variables():
        """Read variable information"""
        names, units, abbrs = [], [], []
        with open(os.path.join(os.path.dirname(__file__),
                  'mascaret_variables_fr.csv'), newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            for row in reader:
                names.append(row['varname'])
                units.append(row['unit'])
                abbrs.append(row['abbr'])
        return names, abbrs, units

    def _read_binary_variables(self):
        variables = []
        while True:
            # First Fortran tag
            self._file.read(4)
            # Read lines to get variables name
            [line] = struct.unpack(self.endians+'4s', self._file.read(4))
            check = line.strip()
            # Ending Fortran tag
            self._file.read(4)
            if len(check) != 0:
                if check.decode('utf-8') == 'FIN':
                    break
                else:
                    variables.append(check.decode('utf-8'))
        return variables

    def _read_binary_header(self):
        """Read header with variables"""

        # Three commentaries block
        for i in range(3):
            # First Fortran tag
            self._file.read(4)
            # Reading commentary (unpack is returning a tuple)
            [line] = struct.unpack(self.endians+'72s', self._file.read(72))
            # If commentaries are not empty, they are printing
            if len(line.split()) != 0:
                self.logger.debug(line)
            # Ending Fortran tag
            self._file.read(4)

        self._read_binary_variables()

        # Number of reaches block
        # First Fortran tag
        self._file.read(4)
        # Number of reaches
        [self.nreaches] = struct.unpack(self.endians+'i', self._file.read(4))
        # Number of reaches repeated
        self._file.read(4)
        # Ending Fortran tag
        self._file.read(4)

        # First and last section of reaches block
        # Defining format to read first and last point (depends on the number
        # of reaches)
        fmt = self.endians+'%ii' % self.nreaches

        # First Fortran tag
        self._file.read(4)
        # Read index of first section (or point) of reaches
        self.reach_first_points = \
            struct.unpack(fmt, self._file.read(self.nreaches*4))
        # Ending Fortran tag
        self._file.read(4)

        # First Fortran tag
        self._file.read(4)
        # Read index of last section (or point) of reaches
        self.reach_last_points = \
            struct.unpack(fmt, self._file.read(self.nreaches*4))
        # Ending Fortran tag
        self._file.read(4)

        # Read variables
        self.list_variables_ind = self._read_binary_variables()

        self._file.read(4)
        [self.nsections] = struct.unpack(self.endians+'i', self._file.read(4))
        self._file.read(4)  # nsections again...
        self._file.read(4)

        self.position_var_ind = self._file.tell()

        self.res_variables_ind = []

        for _ in self.list_variables_ind:
            self._file.read(4)
            self.res_variables_ind.append(
                list(struct.unpack(self.endians+str(self.nsections)+'f',
                     self._file.read(4*self.nsections))))
            self._file.read(4)

        self.list_variables_dep = self._read_binary_variables()

        self.varnames = self.list_variables_dep

        for i in self.list_variables_ind:
            self.varnames.append(i)

        self._position_first_frame = self._file.tell()

        try:
            names, abbrs, units = self._read_dico_variables()
            for i, varname in enumerate(self.varnames):
                var_index = abbrs.index(varname)
                self.add_variable(names[var_index], units[var_index], varname)
        except FileNotFoundError:
            self.logger.warning('Mascaret dico file is missing !')
            for i, varname in enumerate(self.varnames):
                self.add_variable('Rub_long_name_unknown_' + str(i),
                                  'Rub_unit_unknown_' + str(i), varname)

    def get_reaches(self):
        """Read geometry"""
        self._file.seek(self.position_var_ind)

        last_point = 0

        self._file.read(4)
        # Read  time independent variables
        res_variables_ind = \
            list(struct.unpack(self.endians + str(self.nsections) + 'f',
                 self._file.read(4 * self.nsections)))
        self._file.read(4)

        for i in range(self.nreaches):
            # Reaches ID starting at 1
            reach = Reach(i+1)
            for j in range(self.reach_first_points[i],
                           self.reach_last_points[i]+1):
                pk = res_variables_ind[j-1]
                # Dictid header witonary index starting at 1
                reach.add_section(Section(j, pk))

            last_point += j
            # Dictionary index starting at 1
            self._reaches[i+1] = reach

    def get_values(self, record, vars_indexes=None):
        """
        Get values for all variables for a give time index

        @param record (int) time index
        @param vars_indexes (list) List of variable names
        @return (dict) dict of 2D array with reach.id as key
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        res = OrderedDict()
        requested_time = self.times[record]
        self._file.seek(self._times_pos[record])

        values = []
        size = self.nvar - len(self.list_variables_ind)
        for i in range(size):
            self._file.read(4)
            values.append((struct.unpack(self.endians + str(self.nsections) +
                           'f', self._file.read(4*self.nsections))))
            self._file.read(4)
        for i, _ in enumerate(self.list_variables_ind):
            values.append(self.res_variables_ind[i][:])

        all_values = []
        for i in range(self.nsections):
            all_values.append(np.array(values)[:, i][vars_indexes])

        for i, reach in self.reaches.items():
            # shift of index because reach dictionary keys begin at 1 and not 0
            res[reach.id] = \
                np.asarray(all_values[self.reach_first_points[i-1]-1:
                                      self.reach_last_points[i-1]])

        return res

    def get_series(self, reach_id, section_id, vars_indexes=None):
        """
        Get values for all variables for a give reach index and a given section
        index
        @param reach_id (int) reach index
        @param section_id (int) section index
        @param vars_indexes (list) List of variable names
        @return (numpy.array)
        """
        # FIXME: not consistant with Opthyca.get_series =>
        # it should use reach_id argument
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        res = []
        values = []

        var_names = [self.varnames_dict['abbr'][i] for i in vars_indexes]

        for i in var_names:
            if i in self.list_variables_ind:
                var_ind_index = self.list_variables_ind.index(i)

                for j in range(self.ntimestep):
                    values.append(
                        self.res_variables_ind[var_ind_index][section_id])

            else:
                var_dep_index = self.list_variables_dep.index(i)
                for j in range(self.ntimestep):
                    self._file.seek(self._times_pos[j])
                    self._file.read((8+4*self.nsections)*var_dep_index + 4)
                    values.append(np.array(
                        struct.unpack(self.endians + str(self.nsections) + 'f',
                                      self._file.read(4 * self.nsections))
                                           )[section_id]
                                  )
                    self._file.read(4)

            res.append(values)

        return np.array(res).T

    def get_time(self):
        self._move_to_first_frame()
        binary_frame_size = (self.nvar - len(self.list_variables_ind)) *\
            (self.nsections * 4 + 8) + 48
        nb_frames = (self._size_file - self._position_first_frame) //\
            binary_frame_size
        for i in range(nb_frames):
            # skipping 2 integers (frame number x2) + 3 fortran tags
            self._file.read(20)
            [time] = struct.unpack(self.endians+'f', self._file.read(4))
            self._times.append(time)
            self._ntimestep += 1
            self._file.read(24)
            self._times_pos.append(self._file.tell())
            self._file.seek(self._position_first_frame + (i + 1) *
                            binary_frame_size)
        if nb_frames != self._ntimestep:
            self.error("Number of frames is not consistant!")

    def get_endians(self):
        pos_init = self._file.tell()
        self._file.seek(0)

        # first tag
        endian_test = self._file.read(4)

        test_1 = struct.unpack('<i', endian_test)[0]
        self._file.seek(pos_init)
        test_2 = struct.unpack('>i', endian_test)[0]

        if test_1 == 72:
            self._endians = '<'
        elif test_2 == 72:
            self._endians = '>'
        else:
            self.error("Size and alignment of the binary file is neither\
                        little-endian nor big-endian"
                       " or the file is an ASCII File")

        self._file.seek(pos_init)


#  Question : est-il possible de sortir cette classe pour la mettre
#  dans un fichier courlis_file.py?
class ListingCourlis(MascaretFileParent):

    def __init__(self, file_name, access='r', log_lvl='INFO'):
        """
        Constructor for ListingCourlis file

        @param file_name (string) Name of the file
        @param access (string) Access to the file ('r' for read 'w' for write)
        @param log_lvl (string) Logger level

        Attributs specific to ListingCourlis:
        - fformat
        - _size_file
        - _endians
        """
        super().__init__(file_name, access=access + 'b', log_lvl=log_lvl)
        self.fformat = 'listingcourlis'
        self._size_file = os.path.getsize(self.file_name)
        self._endians = ''
        self._sediment_budget_string_res = "Sediment budget not yet extracted"

        self._read_first_time_step()

    @property
    def endians(self):
        """
        Returns the endians of the file
        """
        if not self._endians:
            self.get_endians()
        return self._endians

    @staticmethod
    def _read_dico_variables():
        """Read variable information"""
        names, units, abbrs = [], [], []
        with open(os.path.join(os.path.dirname(__file__),
                  'listing_courlis_variables_fr.csv'), newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            for row in reader:
                names.append(row['varname'])
                units.append(row['unit'])
                abbrs.append(row['abbr'])
        return names, abbrs, units

    def _read_first_time_step(self):
        """
        Reading first time step
        """
        # first timestep pos (first timestep has a different format)
        self._times_pos.append(self._file.tell())

        self._file.read(4)
        res = list(struct.unpack('>17d', self._file.read(17*8)))
        self._file.read(4)

        self._first_timestep = res[0]
        self.nreaches = 1

        end = 0
        count = 1
        prof = 0

        # Search the number of sections and the size of a frame
        try:
            # Initialize the shape of res
            self._file.read(4)
            res = np.vstack((res, list(struct.unpack(self.endians+'17d',
                             self._file.read(17*8)))))
            self._file.read(4)

            while res[count, 1] != 3999:

                file_pos = self._file.tell()
                self._file.read(4)
                res = np.vstack((res, list(struct.unpack(self.endians+'17d',
                                self._file.read(17*8)))))
                self._file.read(4)

                count += 1

                if (res[count, 1] < prof) \
                        and (end == 0):
                    self._res = np.delete(res, -1, 0)
                    self.nsections = count
                    self.reach_first_points = [1]
                    self.reach_last_points = [count]
                    end = 1
                    # first frame corresponds to the second timestep
                    # because the format of the first timestep is
                    # different (no sediment budget)
                    self._position_first_frame = file_pos
                    self._section_idx = self._res[:, 1]
                    self._section_pk = self._res[:, 2]

                if(res[count, 1] == 999):
                    self._position_first_sediment_budget_frame = file_pos
                    self._idx_first_sediment_budget_frame = count

                prof = res[count, 1]

        except(struct.error):
            traceback.print_exc(file=sys.stdout)
            pass

        frame_2_pos_end = self._file.tell()
        self._binary_frame_size = frame_2_pos_end - self._position_first_frame
        self._binary_sediment_budget_frame_size = frame_2_pos_end - \
            self._position_first_sediment_budget_frame
        # It is minus 2 (and not 3 for 999, 1999 and 3999 rows,
        # because Python idx begins to 0
        self.nlayers = count - self._idx_first_sediment_budget_frame - 2

        try:
            names, abbrs, units = self._read_dico_variables()
            for i, varname in enumerate(names):
                self.add_variable(varname, units[i], abbrs[i])

        except FileNotFoundError:
            self.logger.warning('ListingCourlis dico file is missing !')
            for i in range(15):
                self.add_variable('ListingCourlis_long_name_unknown_' + str(i),
                                  'ListingCourlis_unit_unknown_' + str(i),
                                  'ListingCourlis_short_name_unknown_'+str(i))

    def get_reaches(self):
        """Read geometry"""

        # Only one reach in Courlis
        reach = Reach(1)
        for i in range(len(self._section_idx)):
            # Dictionary index starting at 1
            reach.add_section(
                              Section(
                                      self._section_idx[i],
                                      self._section_pk[i]
                                     )
                             )

        self._reaches[1] = reach

    def get_values(self, record, vars_indexes=None):
        """
        Get values for all variables for a give time index

        @param record (int) time index
        @param vars_indexes (list) List of variable names
        @return (dict) dict of 2D array with reach.id as key
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        res = OrderedDict()
        requested_time = self.times[record]
        self._file.seek(self._times_pos[record])

        values = []
        for i in range(self.nsections):
            self._file.read(4)
            values.append((
                struct.unpack(self.endians + '17d', self._file.read(8*17))))
            self._file.read(4)

        values = np.array(values).T
        values = values[3:, :]
        values = np.insert(
                           values,
                           14,
                           values[2, :] * values[3, :],
                           axis=0
                           )
        all_values = []
        for i in range(self.nsections):
            all_values.append(values[:, i][vars_indexes])

        for i, reach in self.reaches.items():
            # shift of index because reach dictionary keys begin at 1 and not 0
            res[reach.id] = \
                np.asarray(all_values[self.reach_first_points[i-1]-1:
                                      self.reach_last_points[i-1]])

        return res

    def get_series(self, reach_id, section_id, vars_indexes=None):
        """
        Get values for all variables for a give reach index and
            a given section index
        @param reach_id (int) reach index
        @param section_id (int) section index
        @param vars_indexes (list) List of variable names
        @return (numpy.array)
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        all_values = []
        section_idx = self.reaches[reach_id].get_section_idx(section_id)
        for record in range(self.ntimestep):
            tmp = self.get_values(record)
            all_values.append(tmp[reach_id][section_idx, vars_indexes])

        return np.array(all_values)

    def get_time(self):
        self._move_to_first_frame()
        nb_frames = (self._size_file - self._position_first_frame) //\
            self._binary_frame_size

        self._times.append(self._first_timestep)
        self._ntimestep += 1

        for i in range(nb_frames):
            self._times_pos.append(self._file.tell())

            self._file.read(4)
            res_first_section = list(struct.unpack(self.endians+'17d',
                                                   self._file.read(17*8)))
            self._file.read(4)

            if res_first_section[0] == self._times[i-1]:
                res_first_section[0] += 0.000001
                self.logger.warning(' Two timesteps, n°{} and n°{}, are the\
                same (often the case for the two first timesteps with sarap)'
                                    ', the second has been augmented\
                                     of +1.e-6s'.format(i, i+1))

            self._times.append(res_first_section[0])
            self._ntimestep += 1
            self._file.seek(self._position_first_frame + (i + 1) *
                            self._binary_frame_size)

    def sediment_budget(self, pandas=False):

        sediment_budget = []

        for i in range(1, self.ntimestep):

            self._file.seek(
                self._position_first_frame +
                i * self._binary_frame_size -
                self._binary_sediment_budget_frame_size
                )

            for j in range(self.nlayers + 3):
                self._file.read(4)
                sediment_budget.append(list(struct.unpack(self.endians+'17d',
                                       self._file.read(17*8))))
                self._file.read(4)

        sediment_budget = np.array(sediment_budget)

        header = ["Timestep",
                  "Sediment mass in suspension",
                  "Initial sediment mass in suspension in the bief",
                  "Sediment mass entering in the bief",
                  "Sediment mass leaving out the bief",
                  "Sediment mass in suspension",
                  "Deposited sediment mass in the bief",
                  ]

        header_total = header.copy()
        header_total.append("Volume sediment variation")

        header.insert(1, "Erosion flux in bief")
        header.insert(2, "Sediment flux entering in the bief")
        header.insert(3, "Sediment flux leaving in the bief")
        for i in range(self.nlayers):
            header.insert(i + 4, "Mass deposited in layer n°"+str(i))
            header_total.insert(i + 1, "Mass deposited in layer n°" + str(i))
        header.insert(self.nlayers + 5, "Relative error for mass budget")
        header.insert(self.nlayers + 11, "Error on mass")
        header.insert(self.nlayers + 12, "Global relative error")

        sand = []
        mud = []
        total = []

        for i in range(self.ntimestep - 1):
            gap = i * (self.nlayers + 3)
            gap_2 = gap + self.nlayers + 1

            mud_temp = [
                        self.times[i+1],
                        sediment_budget[gap, 2],
                        sediment_budget[gap, 4],
                        sediment_budget[gap, 6],
                        sediment_budget[gap_2, 2],
                        sediment_budget[gap_2, 5],
                        sediment_budget[gap_2 + 1, 2],
                        sediment_budget[gap_2 + 1, 4],
                        sediment_budget[gap_2 + 1, 6],
                        sediment_budget[gap_2 + 1, 8],
                        sediment_budget[gap_2 + 1, 10],
                        sediment_budget[gap_2 + 1, 12],
                        sediment_budget[gap_2 + 1, 14]
                       ]

            sand_temp = [
                        self.times[i+1],
                        sediment_budget[gap, 3],
                        sediment_budget[gap, 5],
                        sediment_budget[gap, 7],
                        sediment_budget[gap_2, 3],
                        sediment_budget[gap_2, 6],
                        sediment_budget[gap_2 + 1, 3],
                        sediment_budget[gap_2 + 1, 5],
                        sediment_budget[gap_2 + 1, 7],
                        sediment_budget[gap_2 + 1, 9],
                        sediment_budget[gap_2 + 1, 11],
                        sediment_budget[gap_2 + 1, 13],
                        sediment_budget[gap_2 + 1, 15]
                       ]

            total_temp = [
                        self.times[i+1],
                        sediment_budget[gap_2, 4],
                        sediment_budget[gap_2 + 1, 2] +
                        sediment_budget[gap_2 + 1, 3],
                        sediment_budget[gap_2 + 1, 4] +
                        sediment_budget[gap_2 + 1, 5],
                        sediment_budget[gap_2 + 1, 6] +
                        sediment_budget[gap_2 + 1, 7],
                        sediment_budget[gap_2 + 1, 8] +
                        sediment_budget[gap_2 + 1, 9],
                        sediment_budget[gap_2 + 1, 10] +
                        sediment_budget[gap_2 + 1, 11],
                        sediment_budget[gap_2 + 1, 16]
                       ]

            for j in range(self.nlayers):
                mud_temp.insert(j + 4, sediment_budget[gap + j + 1, 3])
                sand_temp.insert(j + 4, sediment_budget[gap + j + 1, 4])
                total_temp.insert(j + 1, sediment_budget[gap + j + 1, 5])

            mud.append(mud_temp)
            sand.append(sand_temp)
            total.append(total_temp)

        if pandas:
            try:
                import pandas as pd

                mud_df = pd.DataFrame(mud, columns=header)
                sand_df = pd.DataFrame(sand, columns=header)
                total_df = pd.DataFrame(total, columns=header_total)

                return mud_df, sand_df, total_df

            except ImportError:
                print("Pandas module is not available")
                print("Hence, returns objects are numpy array\n")

                # return np.array(mud), np.array(sand), np.array(total)
                return sediment_budget

        else:
            # return np.array(mud), np.array(sand), np.array(total)
            return sediment_budget

    def export_sediment_budget_to_csv(self, outfile_name):
        sediment_budget = self.sediment_budget()
        self._sediment_budget_string_res = ""

        for i in range(0, self.ntimestep - 1):
            gap = i * (self.nlayers + 3)
            gap_2 = gap + self.nlayers + 1
            string_layers = ""

            for j in range(self.nlayers):
                string_layers += \
                    "   Mass deposited in layer n° ;"\
                    + str(j + 1) + 11 * " " + ";" + \
                    "{:16.3f}".format(sediment_budget[gap + j + 1, 3]) +\
                    ";" + " kg" + 10 * " " + ";" + \
                    "{:16.3f}".format(sediment_budget[gap + j + 1, 4]) +\
                    ";" + " kg" + 10 * " " + ";" + \
                    "{:16.3f}".format(sediment_budget[gap + j + 1, 5]) +\
                    ";" + " kg\n"

            self._sediment_budget_string_res += \
                "Timestep : ; {:16.8f}".format(self.times[i + 1]) +\
                ";" + (20 + 9) * " " + \
                "Mud" + ";" + (16 + 10) * " " + "Sand" + (10 + 10) *\
                " " + ";" + "Total\n\n" + \
                "   Erosion flux in the bief            " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap, 2]) + ";" +\
                " kg/s" + 8 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap, 3]) + ";" +\
                " kg/s\n" + \
                "   Sediment flux entering in the bief  " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap, 4]) + ";" +\
                " kg/s" + 8 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap, 5]) + ";" +\
                " kg/s\n" + \
                "   Sediment flux leaving out the bief  " + 4 *\
                " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap, 6]) + ";" +\
                " kg/s" + 8 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap, 7]) + ";"\
                + " kg/s\n\n" \
                " -From t = 0 s to t = {:16.8f} s\n\n"\
                .format(self.times[i + 1]) + \
                string_layers + \
                "   Sediment mass in suspension         " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2, 2]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2, 3]) +\
                ";" + " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2, 4]) +\
                ";" + " kg\n" + \
                "   Relative error for mass budget      " + 4 *\
                " " + ";" + \
                "{:16.8f}".format(sediment_budget[gap_2, 5]) + ";"\
                + 13 * " " + ";" +  \
                "{:16.8f}".format(sediment_budget[gap_2, 6]) + ";\n\n" +\
                " -Global budget mass\n\n" + \
                "   Initial sediment mass in suspension in the bief   " + 4\
                * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 2]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 3]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 2] +
                                  sediment_budget[gap_2 + 1, 3]) +\
                ";" + " kg\n" + \
                "   Sediment mass entering in the bief  " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 4]) +\
                ";" + " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 5]) + ";"\
                + " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 4] +
                                  sediment_budget[gap_2 + 1, 5]) +\
                ";" + " kg\n" + \
                "   Sediment mass leaving out the bief  " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 6]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 7]) + ";"\
                + " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 6] +
                                  sediment_budget[gap_2 + 1, 7]) +\
                ";" + " kg\n" + \
                "   Sediment mass in suspension         " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 8]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 9]) + ";"\
                + " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 8] +
                                  sediment_budget[gap_2 + 1, 9]) + ";"\
                + " kg\n" + \
                "   Deposited sediment mass in the bief " + 4 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 10]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 11]) + ";"\
                + " kg" + 10 * " " + ";" + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 10] +
                                  sediment_budget[gap_2 + 1, 11]) + ";"\
                + " kg\n" + \
                "   Error on mass                       " + 4 * " " + ";" + \
                "{:16.5f}".format(sediment_budget[gap_2 + 1, 12]) + ";" +\
                " kg" + 10 * " " + ";" + \
                "{:16.5f}".format(sediment_budget[gap_2 + 1, 13]) + ";" +\
                " kg" + 10 * " " + ";\n" + \
                "   Global relative error               " + 4 * " " + ";" + \
                "{:16.8f}".format(sediment_budget[gap_2 + 1, 14]) + ";" +\
                13 * " " + ";" +  \
                "{:16.8f}".format(sediment_budget[gap_2 + 1, 15]) + ";\n\n" + \
                " -Volume sediment variation from t = 0 s :\
                    ;{:16.3f}; m3\n\n\n".format(sediment_budget[gap_2 + 1, 16])

        outfile = open(outfile_name, 'w')
        outfile.write(self._sediment_budget_string_res)
        outfile.close()

    def export_sediment_budget_to_txt(self, outfile_name):
        sediment_budget = self.sediment_budget()
        self._sediment_budget_string_res = ""

        for i in range(0, self.ntimestep - 1):
            gap = i * (self.nlayers + 3)
            gap_2 = gap + self.nlayers + 1
            string_layers = ""

            for j in range(self.nlayers):
                string_layers += \
                    "   Mass deposited in layer n° " + str(j + 1) + 12 * " " +\
                    "{:16.3f}".format(sediment_budget[gap + j + 1, 3]) +\
                    " kg" + 10 * " " + \
                    "{:16.3f}".format(sediment_budget[gap + j + 1, 4]) +\
                    " kg" + 10 * " " + \
                    "{:16.3f}".format(sediment_budget[gap + j + 1, 5]) +\
                    " kg\n"

            self._sediment_budget_string_res += \
                "Timestep :  {:16.8f}"\
                .format(self.times[i + 1]) + (20 + 9) * " " + \
                "Mud" + (16 + 10) * " " + "Sand" + (10 + 10) * " " +\
                "Total\n\n" + \
                "   Erosion flux in the bief            " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap, 2])\
                + " kg/s" + 8 * " " + \
                "{:16.3f}".format(sediment_budget[gap, 3]) + " kg/s\n" + \
                "   Sediment flux entering in the bief  " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap, 4]) +\
                " kg/s" + 8 * " " + \
                "{:16.3f}".format(sediment_budget[gap, 5]) + " kg/s\n" + \
                "   Sediment flux leaving out the bief  " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap, 6]) +\
                " kg/s" + 8 * " " + \
                "{:16.3f}".format(sediment_budget[gap, 7]) + " kg/s\n\n" \
                " -From t = 0 s to t = {:16.8f} s\n\n"\
                .format(self.times[i + 1]) + \
                string_layers + \
                "   Sediment mass in suspension         " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2, 2]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2, 3]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2, 4]) +\
                " kg\n" + \
                "   Relative error for mass budget      " + 4 * " " + \
                "{:16.8f}".format(sediment_budget[gap_2, 5]) + 13 * " " +  \
                "{:16.8f}".format(sediment_budget[gap_2, 6]) + "\n\n" +\
                " -Global budget mass\n\n" + \
                "   Initial sediment mass in suspension in the bief   "\
                + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 2])\
                + " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 3])\
                + " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 2]
                                  + sediment_budget[gap_2 + 1, 3])\
                + " kg\n" + \
                "   Sediment mass entering in the bief  " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 4])\
                + " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 5])\
                + " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 4]
                                  + sediment_budget[gap_2 + 1, 5])\
                + " kg\n" + \
                "   Sediment mass leaving out the bief  " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 6]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 7]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 6] +
                                  sediment_budget[gap_2 + 1, 7]) + " kg\n" + \
                "   Sediment mass in suspension         " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 8]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 9]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 8] +
                                  sediment_budget[gap_2 + 1, 9]) + " kg\n" + \
                "   Deposited sediment mass in the bief " + 4 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 10]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 11]) +\
                " kg" + 10 * " " + \
                "{:16.3f}".format(sediment_budget[gap_2 + 1, 10] +
                                  sediment_budget[gap_2 + 1, 11]) + " kg\n" + \
                "   Error on mass                       " + 4 * " " + \
                "{:16.5f}".format(sediment_budget[gap_2 + 1, 12]) + \
                " kg" + 10 * " " + \
                "{:16.5f}".format(sediment_budget[gap_2 + 1, 13]) + \
                " kg" + 10 * " " + "\n" + \
                "   Global relative error               " + 4 * " " + \
                "{:16.8f}".format(sediment_budget[gap_2 + 1, 14]) + \
                13 * " " +  \
                "{:16.8f}".format(sediment_budget[gap_2 + 1, 15]) + "\n\n" + \
                " -Volume sediment variation from t = 0 s : {:16.3f} m3\n\n\n"\
                .format(sediment_budget[gap_2 + 1, 16])

        outfile = open(outfile_name, 'w')
        outfile.write(self._sediment_budget_string_res)
        outfile.close()

    def get_endians(self):
        pos_init = self._file.tell()
        self._file.seek(0)

        # first tag
        # Question Yoann: i ou d pareil?
        endian_test = self._file.read(4)
        test_1 = struct.unpack('<i', endian_test)[0]
        self._file.seek(pos_init)
        test_2 = struct.unpack('>i', endian_test)[0]

        if test_1 == 136:
            self._endians = '<'
        elif test_2 == 136:
            self._endians = '>'
        else:
            self.error("Size and alignment of the binary file is neither\
                        little-endian nor big-endian"
                       " or the file is an ASCII File")

        self._file.seek(pos_init)


class ptravers(MascaretFileParent):

    def __init__(self, file_name, access='r', log_lvl='INFO'):
        """
        Constructor for ptravers Courlis result file

        @param file_name (string) Name of the file
        @param access (string) Access to the file ('r' for read 'w' for write)
        @param log_lvl (string) Logger level

        Attributs specified to ptravers:
        - fformat

        """
        super().__init__(file_name, access=access, log_lvl=log_lvl)
        self.fformat = 'ptravers'
        self._section_pk = []

        self.add_variable('Free surface', 'm', 'Z')
        self.add_variable('Left wet abscissa', 'm', 'LWA')
        self.add_variable('Right wet abscissa', 'm', 'RWA')

        self._read_ptravers_header()

    def read_line(self):
        return self._file.readline().rstrip('\n')

    def _read_ptravers_header(self):
        # read the header of ptravers File

        # skip 2 first lines
        [self.read_line() for i in range(2)]

        # get number of X-proFile
        self.nsections = int(self.read_line().split()[-1:][0])

        # Courlis (and so, ptravers) handle only one reach
        self.reach_first_points = [1]
        self.reach_last_points = [self.nsections + 1]

        # read number of inner points in each X-profile
        self._profile_points = []
        for i in range(self.nsections):
            self._profile_points.append(int(self.read_line().split()[-1:][0]))

        # get number of variables wrote in ptravers
        self._numbersectionvar = int(self.read_line().split()[-1:][0])

        # ptravers file does not contains variable long names
        ptravers_var_dict = {'DXSC': 'Profile abscissa',
                             'ZREF': 'Elevation',
                             'ZDUR': 'Hard bottom elevation',
                             'TauH': 'Local shear stress',
                             'TauE': 'Effective shear stress',
                             'Ceq': 'Equilibrium sand concentration'}

        # read name of variables
        layer_number = 0
        for i in range(self._numbersectionvar):
            line = self.read_line()
            var_abbr = line[:10]
            var_abbr = var_abbr.strip()
            unit = line[10:]
            unit = unit.strip()

            varname = ptravers_var_dict[var_abbr[:4]]
            if var_abbr == 'ZREF':
                if layer_number == 0:
                    var_abbr = 'ZBOT'
                    varname = 'River bottom elevation'
                else:
                    var_abbr += '_' + str(layer_number)
                    varname = 'Layer interface elevation ' + str(layer_number)
                layer_number += 1
            self.add_section_variable(varname, unit, var_abbr)
        self._end_header = self._file.tell()

    def get_reaches(self):
        """Read geometry for ptravers"""

        # Only one reach in Courlis
        reach = Reach(1)
        self._file.seek(self._end_header)

        # skip time in the first X-profile header
        self.read_line()

        for i in range(self.nsections):
            # read X-profile abscissa
            self._section_pk.append(float(self.read_line().split()[-1:][0]))
            # skip free surface and left/right wet abscissa
            self.read_line()
            self.read_line()

            for j in range(self._profile_points[i]):
                self.read_line()

            # Dictionary index starting at 1
            section = Section(i + 1, self._section_pk[i])
            section.nb_points = self._profile_points[i]

            reach.add_section(section)

        self._reaches[1] = reach

    def get_values(self,
                   record,
                   vars_indexes=None,
                   get_section_values=False,
                   section_vars_indexes=None):
        """
        Get values for all variables for a give time index
        and also, for all section variable for a give list of indexes

        @param record (int) time index
        @param vars_indexes (list) List of variable names
        @param get_section_values (boolean) If True get section values as well
        @param section_vars_indexes (list) List of section var index to extract
        @return (dict) dict of 2D array with reach.id as key
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        if get_section_values:
            if section_vars_indexes is None:
                section_vars_indexes = self.section_varnames_dict['id']
            res_section = OrderedDict()
            all_section_values = []

        res = OrderedDict()
        all_values = []

        requested_time = self.times[record]
        self._file.seek(self._times_pos[record])

        # skip time in the first X-profile header
        self.read_line()
        for i in range(self.nsections):

            # skip X-profile abscissa
            self.read_line()

            line = self.read_line()
            free_surface = float(line.split()[-1:][0])

            line = self.read_line()
            left_wet_abscissa = float(line.split()[-2:][0])
            right_wet_abscissa = float(line.split()[-1:][0])

            all_values.append([free_surface, left_wet_abscissa,
                               right_wet_abscissa])

            if get_section_values:
                profile_res = []
                for j in range(self._profile_points[i]):
                    line = self.read_line()
                    profile_res.append(
                                       [float(k) for k in line.split()]
                                      )
                all_section_values.append(np.array(profile_res).T)
            else:
                for j in range(self._profile_points[i]):
                    self.read_line()

        all_values = np.array(all_values).T

        if get_section_values:

            selected_values = []
            selected_section_values = []

            for i in range(self.nsections):
                selected_values.append(all_values[:, i][vars_indexes])
                selected_section_values.append(
                    all_section_values[i][section_vars_indexes])

            for i, reach in self.reaches.items():
                # shift of index because reach dictionary keys
                # begin at 1 and not 0
                res[reach.id] = \
                    np.asarray(selected_values[self.reach_first_points[i - 1]
                               - 1:self.reach_last_points[i - 1]])
                res_section[reach.id] = \
                    selected_section_values[self.reach_first_points[i - 1]
                                            - 1:self.reach_last_points[i - 1]]

            return res, res_section

        else:

            selected_values = []
            for i in range(self.nsections):
                selected_values.append(all_values[:, i][vars_indexes])

            for i, reach in self.reaches.items():
                # shift of index because reach dictionary keys
                # begin at 1 and not 0
                res[reach.id] = \
                    np.asarray(selected_values[self.reach_first_points[i-1]-1:
                                               self.reach_last_points[i-1]])

            return res

    def get_series(self, reach_id, section_id, vars_indexes=None):
        """
        Get values for all variables for a give reach index and
            a given section index
        @param reach_id (int) reach index
        @param section_id (int) section index
        @param vars_indexes (list) List of variable names
        @return (numpy.array)
        """
        if vars_indexes is None:
            vars_indexes = self.varnames_dict['id']

        all_values = []
        section_idx = self.reaches[reach_id].get_section_idx(section_id)
        for record in range(self.ntimestep):
            tmp = self.get_values(record)
            all_values.append(tmp[reach_id][section_idx, vars_indexes])

        return np.array(all_values)

    def get_profile_series(self, reach_id, section_id,
                           vars_section_indexes=None):
        """
        Get values for all section variables for a give reach index
            and a given section index
        @param reach_id (int) reach index
        @param section_id (int) section index
        @param vars_section_indexes (list) List of section variable names
        @return (numpy.array)
        """
        pass

    def get_time(self):
        """
        Initialize time variables
        """

        # if reaches has not be initialized, it is does before seeking
        # in the file
        # otherwise everything is shifted
        reaches = self.reaches.items()
        self._file.seek(self._end_header)

        pos = self._file.tell()
        line = self.read_line()

        while line != '':
            time = float(line.split()[-1:][0])
            self._times.append(time)
            self._times_pos.append(pos)
            self._ntimestep += 1
            for _, reach in reaches:
                for i in range(reach.nsections):
                    for _ in range(reach.sections[i+1].nb_points + 3):
                        self.read_line()

            pos = self._file.tell()
            line = self.read_line()


def MascaretFile(file_name, fformat=None, access='r', log_lvl='INFO'):
    """
    Generic wrapper of mascaret classes will return the class assocaited with
    fformat

    @param file_name (string) Path of the file
    @param fformat (string) File format ('opt', 'listingcourlis', 'ptravers',
    'rub'), optional (detection from extension)
    @param access (string) Access to the file ('r' for read 'w' for write)
    @param log_lvl (string) Logger level

    @returns Class structrure
    """
    # Determine file format from file extension
    if fformat is None:
        fformat = file_name.split(".")[-1]

    if access != 'r':
        raise NotImplementedError('Write access is not supported yet!')
    if fformat == 'opt':
        return Opthyca(file_name, access=access, log_lvl=log_lvl)
    if fformat == 'rub':
        return Rubens(file_name, access=access, log_lvl=log_lvl)
    if fformat == 'listingcourlis':
        return ListingCourlis(file_name, access=access, log_lvl=log_lvl)
    if fformat == 'ptravers':
        return ptravers(file_name, access=access, log_lvl=log_lvl)

    raise FileNotFoundError(
            'The format of the file "{}" is not recognized, '
            'please use "fformat" argument to indicate your file format'
            .format(fformat))


if __name__ == '__main__':
    # Parse every Mascaret Opthyca and Rubens validation files
    # Parse every Courlis Opthyca, Rubens, ListingCourlis
    # and ptravers validation files
    from utils.files import recursive_glob
    try:
        rub_files = recursive_glob(
            os.path.join(os.environ['HOMETEL'], 'examples', 'mascaret'),
            '*.rub')
        opt_files = recursive_glob(
            os.path.join(os.environ['HOMETEL'], 'examples', 'mascaret'),
            '*.opt')
        LC_files = recursive_glob(os.path.join(os.environ['HOMETEL'],
                                  'examples', 'courlis'), '*.listingcourlis')
        PT_files = recursive_glob(os.path.join(os.environ['HOMETEL'],
                                  'examples', 'courlis'), '*.ptravers')

        for file_name in sorted(rub_files + opt_files + LC_files + PT_files):
            if 'sarap.rub' not in file_name:
                masc_file = MascaretFile(file_name)
            else:
                print('Ascii Rubens used for SARAP kernel\
                    output is not yet handled')
                print('So, ', file_name, ' is not tested \n')

            # Display infos about geometry, variables and frames
            print(file_name)
            print(masc_file.summary())

            # Call get_values on first frame to display
            # maximum value for all reach and variables
            # values = masc_file.get_values(0)
            # for reach_id, array in masc_file.reaches.items():
            #     print(np.amax(values[reach_id], axis=0))
    except MascaretException as e:
        print(str(e))
