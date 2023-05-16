"""
MascaretStudy class
===================
"""
import copy
import json
import logging
import os
import shutil
import sys
from collections import OrderedDict
from datetime import datetime
from xml.etree import ElementTree as ETree
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from telapy.api.masc import Mascaret
matplotlib.use('Agg')


def import_study_settings_from_json(setting_file):
    """Loader of a json configuration file"""

    settings_dict = {}
    with open(os.path.join(os.getcwd(),
                           os.path.basename(setting_file)), 'r',
              encoding='utf-8') as file:
        settings = file.read()
        settings_dict = json.loads(
            settings, object_pairs_hook=OrderedDict)

    settings_dict['files']['json'] = os.path.join(os.getcwd(),
                                                  os.path.basename(
                                                  setting_file))
    return settings_dict


def import_study_settings_from_py(setting_file):
    """Loader of a python configuration file"""

    settings_dict = {}
    python_file = os.path.splitext(os.path.basename(setting_file))[0]
    sys.path.append(os.getcwd())
    imp = __import__(python_file)
    settings_dict = imp.settings_dict
    settings_dict['files']['python_settings'] = os.path.join(os.getcwd(),
                                                             python_file +
                                                             ".py")
    return settings_dict


def multi_eval(fun):
    """Decorator to detect space or unique point."""

    def wrapper_fun(self, x_val, *args, **kwargs):
        """Get evaluation from space or point.

        :return: function evaluation(s)
        :rtype: np.array([n_eval], n_feature)
        """
        if not isinstance(x_val, list):
            x_val = [x_val]
            n_eval = 1
            shape_eval = (-1)
        else:
            n_eval = len(x_val)
            shape_eval = (n_eval, -1)

        f_val = [None] * n_eval

        for i, x_i in enumerate(x_val):
            f_val[i] = fun(self, x_i, *args, **kwargs)

        f_val = np.array(f_val).reshape(shape_eval)
        return f_val

    return wrapper_fun


class BoundaryCondition:
    """Boundary Condition"""

    def __init__(self):
        """Constructor.
        """
        self.name = ''
        self.typecode = 0
        self.type_bc = ''
        self.var_bc = ''
        self.idxloi = 0
        self.size_bc = 0
        self.taxis = []
        self.value = []

    def copy(self):
        """Copy a bc"""
        b_c = BoundaryCondition()
        b_c.name = self.name
        b_c.typecode = self.typecode
        b_c.type_bc = self.type_bc
        b_c.var_bc = self.var_bc
        b_c.idxloi = self.idxloi
        b_c.size_bc = self.size_bc
        b_c.taxis = copy.deepcopy(self.taxis)
        b_c.value = copy.deepcopy(self.value)
        return b_c

    def __repr__(self):
        string = '{}: type {}, idxloi {}, size {}'.format(self.name,
                                                          self.type_bc,
                                                          self.idxloi,
                                                          self.size_bc)
        return string


class MascaretStudy():
    """MascaretStudy API."""

    logger = logging.getLogger(__name__)

    def __init__(self, settings_args, xcas=None, lig=None,
                 log_lvl='CRITICAL', iprint=1,
                 working_directory='study',
                 sub_directories=None):
        """Constructor.

        1. Create an instance of the class :meth:`Mascaret`
           with :meth:`Mascaret`.
        2. Creates an instance of Mascaret
           with :meth:`Mascaret.create_mascaret`.
        3. Reads model files from :dict:`settings`
           with :meth:`MascaretStudy.import_model_from_dict`.
        4. Initializes the model with :meth:`MascaretStudy.initialize_model`.

        :param dict or str settings_args: dictionary or file with settings
         (files, perturbations, output)
        :param dict xcas: dictionary of modification in the xcas file.
        :param str log_lvl: log level.
        :param int iprint: save log and results (=1) or not not (=0).
        :param str working_directory: working_directory
        :param str sub_directories: paths associated to the sub-directories
        """
        if sub_directories is None:
            sub_directories = {'figures': 'figures',
                               'output': 'output',
                               'data': 'data'}
        if log_lvl == 'INFO':
            i_log = logging.INFO
        elif log_lvl == 'DEBUG':
            i_log = logging.DEBUG
        else:
            i_log = logging.CRITICAL
        logging.basicConfig(level=i_log)

        # Store the directories
        self.paths = dict(sub_directories)
        self.paths['root'] = os.getcwd()
        self.paths['study'] = os.path.join(self.paths['root'],
                                           working_directory)
        self.paths['figures'] = os.path.join(self.paths['study'],
                                             self.paths['figures'])
        self.paths['output'] = os.path.join(self.paths['study'],
                                            self.paths['output'])
        self.paths['data'] = os.path.join(self.paths['study'],
                                          self.paths['data'])

        def create_directory(directory):
            """  create_directory function """
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.makedirs(directory)

        create_directory(self.paths['study'])
        create_directory(self.paths['figures'])
        create_directory(self.paths['data'])
        create_directory(self.paths['output'])

        # Move to the output directory
        os.chdir(self.paths['output'])

        # Create an instance of the class :meth:`Mascaret`
        self.masc = Mascaret(log_lvl)

        # Create an instance of MASCARET
        self.masc.create_mascaret(iprint)

        # Boundary conditions type dictionary
        self.dico_bc_type = {1: 'Discharge(t)',
                             2: 'Waterlevel(t)',
                             4: 'Rating Z(Q)',
                             5: 'Rating Q(Z)'}
        self.dico_bc_var = {1: 'Discharge',
                            2: 'Level',
                            4: 'Level',
                            5: 'Discharge'}

        # Move to the output directory
        os.chdir(self.paths['root'])

        # Read settings for study
        if isinstance(settings_args, str):
            try:
                settings = import_study_settings_from_py(settings_args)
            except IOError:
                try:
                    settings = import_study_settings_from_json(settings_args)
                except IOError:
                    raise Exception('settings_args', settings_args,
                                    'is not valid config file')
        elif isinstance(settings_args, dict):
            settings = settings_args
        else:
            raise Exception('settings_args', settings_args,
                            'is not valid, should be dict or .py  or .json')

        # Move to the data directory
        os.chdir(self.paths['data'])

        # Read model files and configuration parameters from a dictionary
        self.import_model_from_dict(settings, xcas)

        # Move to the output directory
        os.chdir(self.paths['output'])

        # Initialize model
        self.initialize_model(lig=lig)
        os.chdir(self.paths['root'])

        # Save the model description
        with open(os.path.join(self.paths['study'], 'model_description.txt'),
                  'w') as file:
            file.write(repr(self))

    def __del__(self):
        """Delete an instance of MascaretStudy."""
        self.logger.info('Deleting instance #{}...'.format(
            self.masc.id_masc.value))
        self.clean_dump()
        del self.masc
        os.chdir(self.paths['root'])

    def clean_dump(self):
        """Free saved states"""
        self.masc.free_all_saved_states()
        self.dump_index = []

    def import_model_from_dict(self, settings_dict, xcas=None):
        """Read model files (.xcas, .geo, .lig, .loi, .dtd)
        from :dict:`settings_dict` import the model and
        store information concerning the spatio-temporal discretization. Read
        also optional parameters such as initial conditions, output,
        boundary conditions, friction coefficient, bathymetry, ...

        Uses :meth:`Mascaret.import_model`.

        :param dict settings_dict: the settings dictionary
        :param dict xcas: dictionary of modification in the xcas file.
        """
        self.settings = settings_dict
        # Convert all values of settings from str to byte
        file_types = []
        file_names = []
        ref_xcas = ""
        for key_val in self.settings['files'].items():
            # try:
            if not isinstance(key_val[1], list):
                value = key_val[1]
                file_names.append(value)
                self.settings['files'][key_val[0]] = os.path.basename(value)
                file_types.append(key_val[0])
                if file_types[-1] not in ['res', 'listing']:
                    os.symlink(os.path.join(self.paths['root'],
                               file_names[-1]),
                               os.path.join(self.paths['data'],
                                            os.path.basename(file_names[-1])))
                    file_names[-1] = os.path.basename(file_names[-1])
                else:
                    self.settings['files'][key_val[0]] = \
                                os.path.join(self.paths['output'], value)
                    file_names[-1] = os.path.join(self.paths['output'],
                                                  file_names[-1])

                if (file_types[-1] == 'xcas') and (xcas is not None):
                    ref_xcas = os.path.basename(key_val[1])
                    # file_names[-1] = 'modified_'+key_val[1]
                    file_names[-1] = 'modified_' + \
                                     str(xcas['nbSectionZone']) + \
                                     '_' + ref_xcas
                    path_tmp = os.path.join(self.paths['data'],
                                            'modified_' +
                                            self.settings['files']['xcas'])
                    txt_tmp = self.settings['files']['xcas'] + ' --> ' + \
                        path_tmp

                    self.settings['files']['xcas'] = txt_tmp
                    file_names[-1] = file_names[-1]
            # except AttributeError:  # In case of a list, loop over it
            # except TypeError:
            else:
                for i, sub in enumerate(key_val[1]):
                    sub_value = sub
                    file_names.append(sub_value)
                    self.settings['files'][key_val[0]][i] = \
                        os.path.basename(sub_value)
                    file_types.append(key_val[0])
                    os.symlink(os.path.join(self.paths['root'],
                               file_names[-1]),
                               os.path.join(self.paths['data'],
                               os.path.basename(file_names[-1])))
                    file_names[-1] = os.path.basename(file_names[-1])
                    print('InExcept', file_names[-1])

        # Modify the xcas file using a dictionary.
        if xcas is not None:
            tree = ETree.parse(os.path.join(self.paths['root'], ref_xcas))
            for x_val in xcas:
                tree.find('.//' + x_val).text = str(xcas[x_val])
            with open(os.path.join(self.paths['root'], ref_xcas)) as file:
                content = file.readlines()
            line0 = content[0]
            line1 = content[1]
            path_tmp = 'modified_' + str(xcas['nbSectionZone']) + '_' + \
                       ref_xcas
            with open(os.path.join(self.paths['data'], path_tmp), 'wb') as out:
                #  with open(os.path.join(self.paths['data'],
                # 'modified_'+ref_xcas), 'wb')
                # as out:
                out.write(line0)
                out.write(line1)
                tree.write(out, xml_declaration=False, encoding='ISO-8859-1')

        # Import a model
        self.masc.import_model(file_names, file_types)

        # Store information concerning spatial and temporal discretization
        self.store_simulation_times()
        self.store_model_size()
        self.store_curvilinear_abscissa()
        self.store_cross_sections()
        self.store_xcoord()

        # Build boundary conditions list
        self.boundary_conditions()

        # Change some values (boundary conditions, Strickler coefficients,
        # bathymetry, ...)
        if 'boundary_conditions' in self.settings:
            boundary_conditions = self.settings['boundary_conditions']
            for b_c in boundary_conditions:
                self.change_boundary_condition(b_c)
        if 'friction_coefficients' in self.settings:
            friction_coefficients = self.settings['friction_coefficients']
            for k_s in friction_coefficients:
                if k_s['type'] == 'zone':
                    self.set_zone_friction_minor(k_s)
                else:
                    self.set_friction_minor(k_s)
        if 'friction_coefficients_major' in self.settings:
            friction_coefficients_major = \
                self.settings['friction_coefficients_major']
            for k_s in friction_coefficients_major:
                if k_s['type'] == 'zone':
                    self.set_zone_friction_major(k_s)
                else:
                    self.set_friction_major(k_s)
        if 'bathymetry' in self.settings:
            self.perturbate_cross_sections(self.settings['bathymetry'])

    def store_model_size(self):
        """Store the model size (number of nodes)

        Uses :meth:`Mascaret.get_var_size`.
        """
        self.model_size, _, _ = self.masc.get_var_size('Model.X')

    def store_curvilinear_abscissa(self):
        """Store the curvilinear abscissa of the different nodes.

        Use :meth:`Mascaret.get_var_size` and :meth:`Mascaret.get`.
        """
        self.logger.debug('Getting the value of all Model.X...')
        self.curvilinear_abscissa = [self.masc.get('Model.X', i)
                                     for i in range(self.model_size)]
        self.logger.debug('All Model.X get.')

    def store_simulation_times(self):
        """Get the simulation times (initial time, final time and time step).

        Uses Mascaret Api :meth:`Mascaret.get`.
        """
        self.time_step = self.masc.get('Model.DT')
        self.initial_time = self.masc.get('Model.InitTime')
        self.final_time = self.masc.get('Model.MaxCompTime')

    def store_cross_sections(self):
        """Store the different cross sections.

        Uses :meth:`Mascaret.get_var_size` and :meth:`Mascaret.get`.
        """
        # Model.CrossSection.RelAbs curvilinear absissa for cross sections
        # Model.CrossSection.X y-coordinate on the cross sections
        # nb_secs: number of cross sections (scalar)
        # np_pts: number of points per cross sections (1D array, size nb_secs)
        # bot_x: curvilinear abscissa for cross sections
        # (1D array, size nb_secs)
        # bot_z: bottom level for cross sections (1D array, size nb_secs)
        nb_secs, _, _ = self.masc.get_var_size('Model.CrossSection.RelAbs')
        nb_pts = []
        for sec in range(nb_secs):
            _, idx_nb_pts, _ = self.masc.get_var_size('Model.CrossSection.X',
                                                      sec)
            nb_pts.append(idx_nb_pts)
        bot_z = [self.masc.get('Model.CrossSection.Zbot', sec)
                 for sec in range(nb_secs)]
        bot_s = [self.masc.get('Model.CrossSection.RelAbs', sec, 0)
                 for sec in range(nb_secs)]
        s_val = np.repeat(bot_s, nb_pts, axis=0)
        # Model.CrossSection.Y z-coordinate on the cross sections
        # z: z-coordinate for points on the cross section
        # list_z: z-coordinate for points on all the cross sections
        #         list of lists by section)
        list_z = []
        for sec in range(nb_secs):
            list_sec_z = []
            for idpt in range(nb_pts[sec]):
                z_val = self.masc.get('Model.CrossSection.Y', sec, idpt)
                list_sec_z.append(z_val)
            list_z.append(list_sec_z)

        info = ("Getter Cross Section summary:\n"
                "table X: {}\n"
                "table X bot: {}\n"
                "table Z bot: {}").format(s_val, bot_s, bot_z)

        self.cross_sections = {'bottom': {'abscissa': bot_s, 'level': bot_z},
                               'all': {'abscissa': s_val, 'level': list_z}}

        self.logger.info(info)

    def initialize_model(self, z_init=None, q_init=None, lig=None):
        """Initialize the model:
            * from fields
              with z_init and q_init
            * from constant values ``init_cst`` in :attr:`settings` along
              with ``Q_cst`` and ``Z_cst`` values
            * or from :file:`file.lig` in :attr:`settings`.
        Uses :meth:`Mascaret.init_hydro`
        or :meth:`Mascaret.init_hydro_from_file`.
        """
        # Initialize Mascaret Model from fields
        if q_init is not None and z_init is not None:
            self.masc.init_hydro(z_init, q_init)
        elif 'initial_conditions' in self.settings:
            # Initialize Mascaret Model from constant values
            # (from settings dict)
            q_val = [self.settings['initial_conditions']['Q_cst']] * \
                self.model_size
            z_val = [self.settings['initial_conditions']['Z_cst']] * \
                self.model_size
            self.masc.init_hydro(z_val, q_val)
        elif lig is not None:
            lig = os.path.join(self.paths['root'], lig)
            # Initialize Mascaret from run-time provided file
            self.masc.init_hydro_from_file(lig)
        else:
            # Initialize Mascaret Model from file
            self.masc.init_hydro_from_file(
                os.path.join(self.paths['data'],
                             self.settings['files']['lig']))

    def __repr__(self):
        """Class informations based on settings."""

        string = "======================\n"
        string += "= MODEL DESCRIPTION  =\n"
        string += "======================\n"

        string += "\nDIRECTORIES:\n"

        string += " -- study: {}\n"
        values = [self.paths['study']]

        string += " -- Figures: {}\n"
        values.append(self.paths['figures'])

        string += " -- Data: {}\n"
        values.append(self.paths['data'])

        string += " -- Output: {}\n"
        values.append(self.paths['output'])

        string += "\nMODEL FILES:\n"

        if 'json' in self.settings['files']:
            string += " -- json: {}\n"
            values.append(self.settings['files']['json'])

        if 'python_settings' in self.settings['files']:
            string += " -- py settings: {}\n"
            values.append(self.settings['files']['python_settings'])

        string += " -- xcas: {}\n"
        values.append(self.settings['files']['xcas'])

        string += " -- geo: {}\n"
        values.append(self.settings['files']['geo'])

        if 'casier' in self.settings['files']:
            string += " -- casier: {}\n"
            values.append(self.settings['files']['casier'])

        string += " -- res: {}\n"
        values.append(self.settings['files']['res'])

        if 'casier' in self.settings['files']:
            string += " -- res_casier: {}\n"
            values.append(self.settings['files']['listing_casier'])
            string += " -- res_liaison: {}\n"
            values.append(self.settings['files']['listing_liaison'])

        string += " -- listing: {}\n"
        values.append(self.settings['files']['listing'])

        string += " -- damocle: {}\n"
        values.append(self.settings['files']['damocle'])

        string += " -- lig: {}\n"
        values.append(self.settings['files']['lig'])

        string += " -- loi:\n"
        for f_loi in self.settings['files']['loi']:
            string += '         {}\n'
            values.append(f_loi)

        string += '\nUSER SETTINGS:\n'
        if 'initial_conditions' in self.settings:
            string += " -- Initialize model with constant values:\n"
            string += "       > Q_cst: {}\n"
            values.append(self.settings['initial_conditions']['Q_cst'])
            string += "       > Z_cst: {}\n"
            values.append(self.settings['initial_conditions']['Z_cst'])

        if 'boundary_conditions' in self.settings:
            boundary_conditions = self.settings['boundary_conditions']
            for b_c in boundary_conditions:
                string += " -- Change BC:\n"
                string += "       > Name: {}\n"
                values.append(b_c['name'])
                if 'value' in b_c:
                    string += "       > Value: {}\n"
                    values.append(b_c['value'])
                elif 'multcoeff' in b_c:
                    string += "       > Mult coeff: {}\n"
                    values.append(b_c['multcoeff'])
                elif 'addperturb' in b_c:
                    string += "       > Add perturb: {}\n"
                    values.append(b_c['addperturb'])

        if 'friction_coefficients' in self.settings:
            string += " -- Change the friction coefficient:\n"
            friction_coefficients = self.settings['friction_coefficients']
            for k_s in friction_coefficients:
                if k_s['type'] == 'node':
                    string += "       > Node index: {}\n"
                    values.append(k_s['index'])
                    string += "       > Value: {}\n"
                    values.append(k_s['value'])
                elif k_s['type'] == 'zone':
                    string += "       > Zone index: {}\n"
                    values.append(k_s['index'])
                    string += "       > Value: {}\n"
                    values.append(k_s['value'])

        if 'friction_coefficients_major' in self.settings:
            string += " -- Change the friction coefficient for major:\n"
            friction_coefficients_major = \
                self.settings['friction_coefficients_major']
            for k_s in friction_coefficients_major:
                if k_s['type'] == 'node':
                    string += "       > Node index: {}\n"
                    values.append(k_s['index'])
                    string += "       > Value: {}\n"
                    values.append(k_s['value'])
                elif k_s['type'] == 'zone':
                    string += "       > Zone index: {}\n"
                    values.append(k_s['index'])
                    string += "       > Value: {}\n"
                    values.append(k_s['value'])

        if 'bathymetry' in self.settings:
            string += " -- Change bathymetry:\n"
            if 'section' in self.settings['bathymetry']:
                string += "Cross section index: {}\n"
                values.append(self.settings['bathymetry']['section'])
            else:
                string += "Cross section index: all bathymetry"
            if 'variation' in self.settings['bathymetry']:
                string += "     * Shift: {}\n"
                values.append(self.settings['bathymetry']['variation'])
            else:
                string += "     * Gp covariance: {}\n"
                values.append(self.settings['bathymetry']['covariance'])

        if 'output' in self.settings:
            string += " -- Return the hydraulic state at:\n"
            if 'curv_abs' in self.settings['output']:
                string += "     * Curvilinear abscissa: {}\n"
                values.append(self.settings['output']['curv_abs'])
            elif 'index' in self.settings['output']:
                string += "     * Node index: {}\n"
                values.append(self.settings['output']['node'])
            else:
                string += "     * All nodes\n"

        string += "\nMODEL FEATURES:\n"
        string += " -- Spatial features:\n"
        string += "     * First curvilinear abscissa: {}\n"
        values.append(self.curvilinear_abscissa[0])
        string += "     * Last curvilinear abscissa: {}\n"
        values.append(self.curvilinear_abscissa[-1])
        string += "     * Number of nodes: {}\n"
        values.append(self.model_size)
        string += " -- Temporal features:\n"
        string += "     * Initial time: {}\n"
        values.append(self.initial_time)
        string += "     * Final time: {}\n"
        values.append(self.final_time)
        string += "     * Time step: {}\n"
        values.append(self.time_step)

        return string.format(*values)

    def __call__(self, x_val=None, flag='all', dump_state=None, tstart=None,
                 tend=None):
        """Run the Mascaret study using the user configuration
        (:attr:`settings`) and possibly new values of some parameters such
        as the Strickler coefficient Ks, the boundary conditions Q, ...

        :param dict x_val: input parameters {nameX1: valX1, nameX2: valX2, ...}
        :param list dump_state: list of times in seconds for state dumping
        (reload saved state i with self.masc.set_state(self.dump_index[i])
        :param dump_state: list of time for saving state
        :param tstart: initial time
        :param tend: final time
        :return [list(]dict[)] hydraulic_state: dictionary containing
        curvilinear abscissa, water level and flowrate
        """
        # Pre-treatment: change some parameter values
        if x_val is not None:
            if 'friction_coefficients' in x_val:
                for friction_coefficient in x_val['friction_coefficients']:
                    value = friction_coefficient['value']
                    if friction_coefficient['type'] == 'zone':
                        index = friction_coefficient['index']
                        self.set_zone_friction_minor({'index': index,
                                                      'value': value})
                    else:
                        index = friction_coefficient['index']
                        self.set_friction_minor({'index': index,
                                                 'value': value})
            if 'friction_coefficients_major' in x_val:
                for friction_coefficient_major in \
                        x_val['friction_coefficients_major']:
                    value = friction_coefficient_major['value']
                    if friction_coefficient_major['type'] == 'zone':
                        index = friction_coefficient_major['index']
                        self.set_zone_friction_major({'index': index,
                                                      'value': value})
                    else:
                        index = friction_coefficient_major['index']
                        self.set_friction_major({'index': index,
                                                 'value': value})
            if 'boundary_conditions' in x_val:
                boundary_conditions = x_val['boundary_conditions']
                for boundary_condition in boundary_conditions:
                    self.change_boundary_condition(boundary_condition)

            if 'bathymetry' in x_val:
                print('bathymetry', x_val['bathymetry'])
                self.perturbate_cross_sections(x_val['bathymetry'])

        # Move to the output directory
        os.chdir(self.paths['output'])
        # Treatment: execute Mascaret
        if tstart is not None:
            self.initial_time = tstart
        if tend is not None:
            self.final_time = tend
        self.dump_index = []
        if dump_state is None:
            self.masc.compute(self.initial_time, self.final_time,
                              self.time_step)
            self.masc.error_message()
        else:
            t_0 = self.initial_time
            for ti in dump_state:
                if ti > t_0:
                    self.masc.compute(t_0, ti, self.time_step)
                self.dump_index.append(self.masc.save_state())
                t_0 = ti
            if dump_state[-1] < self.final_time:
                self.masc.compute(dump_state[-1], self.final_time,
                                  self.time_step)

        # Post-treatment: design the output (local or global hydraulic state)
        res = None
        if flag is not None and flag != 'none':

            if 'output' in self.settings:
                if 'curv_abs' in self.settings['output']:
                    hydraulic_state = self.local_hydraulic_state
                elif 'node' in self.settings['output']:
                    hydraulic_state = self.local_hydraulic_state
                else:
                    hydraulic_state = self.global_hydraulic_state
            else:
                hydraulic_state = self.global_hydraulic_state

            if flag == 'all':
                res = hydraulic_state
            elif flag == 'z':
                res = hydraulic_state['z']
            elif flag == 'q':
                res = hydraulic_state['q']
            elif flag == 'h':
                res = hydraulic_state['h']

        # Back to root
        os.chdir(self.paths['root'])
        return res

    @property
    def global_hydraulic_state(self):
        """Get the hydraulic state at each node of the mesh.

        :return: curvilinear abscissa, height, flow rate
        :rtype: dict(list(float), list(float), list(float))
        """
        return {'s': self.curvilinear_abscissa,
                'z': self.global_elevation,
                'q': self.global_discharge,
                'zb': self.global_bathymetry,
                'h': self.global_height}

    @property
    def global_height(self):
        """Get the water height at each node of the mesh.

        :return: Water height at each node of the mesh.
        :rtype: list of floats
        """
        elevation = self.global_elevation
        bathymetry = self.global_bathymetry
        height = [elevation[i] - bathymetry[i] for i in range(self.model_size)]
        return height

    @property
    def global_elevation(self):
        """Get the water elevation at each node of the mesh.

        Use Mascaret Api :meth:`Mascaret.get`.

        :return: Water elevation at each node of the mesh.
        :rtype: list of floats
        """
        self.logger.debug('Getting the value of all State.Z...')
        elevation = [self.masc.get('State.Z', i)
                     for i in range(self.model_size)]
        self.logger.debug('All State.Z get.')
        return elevation

    @property
    def global_discharge(self):
        """Get the discharge at each node of the mesh.

        Use Mascaret Api :meth:`Mascaret.get`.

        :return: Discharge at each node of the mesh.
        :rtype: list of floats
        """
        self.logger.debug('Getting the value of all State.Q...')
        discharge = [self.masc.get('State.Q', i)
                     for i in range(self.model_size)]
        self.logger.debug('All State.Q get.')
        return discharge

    @property
    def global_bathymetry(self):
        """Get the bathymetry at each node of the mesh.

        Use Mascaret Api :meth:`Mascaret.get`.

        :return: Bathymetry at each node of the mesh.
        :rtype: list of floats
        """
        self.logger.debug('Getting the value of all Model.Zbot...')
        bathymetry = [self.masc.get('Model.Zbot', i)
                      for i in range(self.model_size)]
        self.logger.debug('All Model.Zbot get.')
        return bathymetry

    @property
    def local_hydraulic_state(self):
        """Get the hydraulic state at a specific curvilinear abscissa or node.

        :return: curvilinear abscissa, height, flow rate
        :return: dict(float, float, float)
        """
        if 'curv_abs' in self.settings['output']:

            def mean_value(x_val, y_1, y_2, x_1, x_2):
                y_mean = (y_2 * (x_val - x_1) + y_1 * (x_2 - x_val)) / \
                 (x_2 - x_1)
                return y_mean

            curv_abs = self.settings['output']['curv_abs']
            temp = [s > curv_abs for s in self.curvilinear_abscissa]
            id_tick = next(index
                           for index, value in enumerate(temp) if value) - 1
            curv_abs_1 = self.curvilinear_abscissa[id_tick]
            elevation_1 = self.local_elevation(id_tick)
            discharge_1 = self.local_discharge(id_tick)
            bathymetry_1 = self.local_bathymetry(id_tick)
            curv_abs_2 = self.curvilinear_abscissa[id_tick + 1]
            elevation_2 = self.local_elevation(id_tick + 1)
            discharge_2 = self.local_discharge(id_tick + 1)
            bathymetry_2 = self.local_bathymetry(id_tick + 1)
            elevation = mean_value(curv_abs,
                                   elevation_1, elevation_2,
                                   curv_abs_1, curv_abs_2)
            discharge = mean_value(curv_abs,
                                   discharge_1, discharge_2,
                                   curv_abs_1, curv_abs_2)
            bathymetry = mean_value(curv_abs,
                                    bathymetry_1, bathymetry_2,
                                    curv_abs_1, curv_abs_2)
            height = elevation - bathymetry
        else:
            index = self.settings['output']['node']
            curv_abs = self.curvilinear_abscissa[index]
            elevation = self.local_elevation(index)
            discharge = self.local_discharge(index)
            bathymetry = self.local_bathymetry(index)
            height = elevation - bathymetry

        return {'s': curv_abs,
                'z': elevation,
                'q': discharge,
                'zb': bathymetry,
                'h': height}

    def local_interpolation(self, listcurv):
        """Provides interpolation information.

        :return: coordinate, with interpolation indices and weight
        :return: dict(float, float, float)
        """

        def weight(x_val, x_1, x_2):
            if x_1 == x_2:
                w_1 = 0.5
                w_2 = 0.5
            else:
                w_1 = (x_2 - x_val) / (x_2 - x_1)
                w_2 = 1. - w_1
            return w_1, w_2

        output = {}
        for curv_abs in listcurv:
            if curv_abs <= self.curvilinear_abscissa[0]:
                id_tick = 0
                id_tickp1 = 0
            elif curv_abs >= self.curvilinear_abscissa[-1]:
                id_tick = self.model_size - 1
                id_tickp1 = self.model_size - 1
            else:
                temp = [s > curv_abs for s in self.curvilinear_abscissa]
                id_tickp1 = next(index
                                 for index, value in enumerate(temp) if value)
                id_tick = id_tickp1 - 1
            curv_abs_1 = self.curvilinear_abscissa[id_tick]
            curv_abs_2 = self.curvilinear_abscissa[id_tickp1]
            weight_1, weight_2 = weight(curv_abs, curv_abs_1, curv_abs_2)
            output[curv_abs] = {'ind1': (id_tick, weight_1),
                                'ind2': (id_tickp1, weight_2)}
        return output

    def local_discharge(self, index):
        """Get the discharge at a given node.

        Uses :meth:`Mascaret.get`.

        :param float index: Node index.
        :return: discharge at this node.
        :rtype: float.
        """
        self.logger.debug('Getting the value of State.Q...')
        discharge = self.masc.get('State.Q', index)
        self.logger.debug('State.Q get.')
        return discharge

    def local_elevation(self, index):
        """Get the water elevation at a given node.

        Uses :meth:`Mascaret.get`.

        :param float index: Node index.
        :return: Water elevation at this node.
        :rtype: float.
        """
        self.logger.debug('Getting the value of State.Z...')
        self.logger.debug('State.Z get.')
        elevation = self.masc.get('State.Z', index)
        return elevation

    def local_bathymetry(self, index):
        """Get the bathymetry at a given node.

        Uses :meth:`Mascaret.get`.

        :param float index: Node index.
        :return: Bathymetry at this node.
        :rtype: float.
        """
        self.logger.debug('Getting the value of Model.Zbot...')
        self.logger.debug('Model.Zbot get.')
        bathymetry = self.masc.get('Model.Zbot', index)
        return bathymetry

    def boundary_conditions(self):
        """Get the numbers, the names and the indices of all laws

        Uses :meth:`masc.get`
        and :meth:`masc.get_var_size`.

        :return: number of BC, names of BC, indices of BC
        :rtype: int, list(str), list(int)
        """
        number, _, _ = self.masc.get_var_size('Model.Graph.Name')
        # Create a boundary condition dictionary
        self.bc = {}

        for k in range(number):
            my_bc = BoundaryCondition()
            my_bc.idxloi = k
            # AP my_bc.idxloi = indices_names[0][k]-1
            my_bc.typecode = self.masc.get('Model.Graph.Type', my_bc.idxloi)
            if my_bc.typecode not in self.dico_bc_type:
                continue
            my_bc.name = self.masc.get('Model.Graph.Name', my_bc.idxloi)
            my_bc.type_bc = self.dico_bc_type[my_bc.typecode]
            # self.masc.get('Model.Graph.Type', my_bc.idxloi)]
            my_bc.var_bc = self.dico_bc_var[my_bc.typecode]
            # self.masc.get('Model.Graph.Type', my_bc.idxloi)]

            _, my_bc.size_bc, _ = self.masc.get_var_size('Model.Graph.'
                                                         + my_bc.var_bc,
                                                         my_bc.idxloi)
            my_bc.value = np.ones(my_bc.size_bc, float)
            if my_bc.typecode == 1 or my_bc.typecode == 2:
                my_bc.taxis = np.ones(my_bc.size_bc, float)
            for t in range(my_bc.size_bc):
                my_bc.value[t] = self.masc.get('Model.Graph.' +
                                               my_bc.var_bc, my_bc.idxloi, t)
                if my_bc.typecode == 1 or my_bc.typecode == 2:
                    my_bc.taxis[t] = self.masc.get('Model.Graph.Time',
                                                   my_bc.idxloi, t)

            self.bc[my_bc.name] = my_bc
            self.bc_keep = copy.deepcopy(self.bc)

        # self.logger.info(' AVAILABLE CONDITIONS FOR PERTURBATIONS')
        # for x in self.BC:
        #    self.logger.info(self.BC[x].name,' of type ',self.BC[x].typeBC,
        # ' as law ',self.BC[x].idxloi)
        # info = ("AVAILABLE CONDITIONS FOR PERTURBATIONS\n"
        #        "{} of type {} as law {}\n").format(self.BC[].name,
        # self.BC[].typeBC,self.BC[].idxloi)
        # self.logger.info(info)

    def change_boundary_condition(self, b_c, keep=True):
        """Set boundary condition Q(t) or Z(t).
        Use :meth:`Mascaret.get_var_size` and :meth:`Mascaret.set`.
        :param dict b_c: Boundary Condition Q(t) or Z(t)
         ``{'name','value' or 'multcoeff' or 'addperturb' or
          'shift_chronicle'}``
        :param bool keep: 'multcoeff', 'addperturb' and 'shift_chronicle'
         are applied on bc_keep values (True) or bc values (False)
        """
        self.logger.debug('change_boundary_condition')

        try:
            if keep:
                self.bc[b_c['name']] = self.bc_keep[b_c['name']].copy()
            x_val = self.bc[b_c['name']]
        except Exception:
            self.logger.exception(' change_boundary_condition: ' + b_c['name']
                                  + ' not in bound cond list')
            return

        try:
            for t in range(x_val.size_bc):
                # Imposed values perturbations take priority
                if 'value' in b_c:
                    b_c['value'] = np.atleast_1d(b_c['value'])
                    x_val.value[t] = b_c['value'][min(t,
                                                      len(b_c['value']) - 1)]
                # Followed by multiplicative coefficients
                if 'multcoeff' in b_c:
                    b_c['multcoeff'] = np.atleast_1d(b_c['multcoeff'])
                    x_val.value[t] *= b_c['multcoeff'][min(t,
                                                       len(b_c['multcoeff'])
                                                           - 1)]
                # And finally by additive perturbations
                if 'addperturb' in b_c:
                    x_val.value[t] += b_c['addperturb'][min(t,
                                                        len(b_c['addperturb'])
                                                            - 1)]

            if 'shift_chronicle' in b_c:
                if x_val.typecode == 1 or x_val.typecode == 2:
                    x_shift = copy.deepcopy(x_val.value)
                    for t in range(x_val.size_bc):
                        t_shift = x_val.taxis[t] + b_c['shift_chronicle']
                        if b_c['shift_chronicle'] < 0.:
                            if t_shift < 0:
                                continue
                            ind = [k for k, bck in enumerate(x_val.taxis) if
                                   bck < t_shift][-1]
                            if t == x_val.size_bc - 1:
                                x_shift[ind + 1:] = x_val.value[-1]
                            slope = (x_val.value[t] - x_val.value[t - 1]) / \
                                    (1. * x_val.taxis[ind + 1] - 1. *
                                     x_val.taxis[ind])
                            x_shift[ind] = \
                                x_val.value[t - 1] + slope * \
                                (x_val.taxis[ind + 1] - t_shift)
                        else:
                            if t_shift > x_val.taxis[-1]:
                                break
                            ind = [k for k, bck in enumerate(x_val.taxis)
                                   if bck < t_shift][-1]
                            if t == 0:
                                x_shift[:ind + 1] = x_val.value[0]
                            slope = (x_val.value[t + 1] - x_val.value[t]) / \
                                    (1. * x_val.taxis[t + 1] - 1. *
                                     x_val.taxis[t])
                            x_shift[ind + 1] = \
                                x_val.value[t] + slope * (x_val.taxis[ind + 1]
                                                          - t_shift)
                    x_val.value = x_shift
                else:
                    self.logger.warning('Shift chronicle only\
                     applies on Q(t) or Z(t)')

            for t in range(x_val.size_bc):
                self.masc.set('Model.Graph.' + x_val.var_bc, x_val.value[t],
                              x_val.idxloi, t)

            self.logger.debug(' BC ' + x_val.name + ' ' + x_val.type_bc +
                              ' set to {}'.format(x_val.value))

        except Exception:
            self.logger.exception(' change_boundary_condition: ' + b_c['name']
                                  + ' went wrong')

    def get_friction_zone(self, index):
        """Get indices of the beginning and end of a friction zone.

        Uses :meth:`Mascaret.get_var_size` and :meth:`Mascaret.get_int`.

        :param: int index: index of a specific zone

        :return: Indices of beginning and end of a friction zone
        :rtype: int, int
        """
        # zone_start and zone_end start from 1,
        # they are shifted here for python
        zone_start = self.masc.get('Model.FrictionZone.FirstNode', index) - 1
        zone_end = self.masc.get('Model.FrictionZone.LastNode', index) - 1
        return zone_start, zone_end

    def get_zone_friction_minor(self, index):
        """Get minor friction coefficient at a given zone.

        Uses :meth:`get_friction_zone` and :meth:`get_friction_minor`.

        :param int index: zone index
        :return: Friction coefficient at zone
        :rtype: list(float)
        """
        zone_start, zone_end = self.get_friction_zone(index)
        zone_friction = self.get_friction_minor()
        zone_friction = zone_friction[zone_start:zone_end]
        self.logger.debug('Zone Ks minor get.')
        return zone_friction

    def get_zone_friction_major(self, index):
        """Get major friction coefficient at a given zone.

        Uses :meth:`get_friction_zone` and :meth:`get_friction_major`.

        :param int index: zone index
        :return: Friction coefficient at zone
        :rtype: list(float)
        """
        zone_start, zone_end = self.get_friction_zone(index)
        zone_friction = self.get_friction_major()
        zone_friction = zone_friction[zone_start:zone_end]
        self.logger.debug('Zone Ks major get.')
        return zone_friction

    def set_zone_friction_minor(self, ks):
        """Change minor friction coefficient at zone.

        Uses :meth:`get_friction_zone` and :meth:`set_friction_minor`.
        Does not use get_zone_friction_minor but loop on set_friction_minor
        as perturbation to ks is not additive.

        :param dict ks: Friction coefficient at zone ``{'index','value'}``
        """
        zone_index, zone_value = ks['index'], ks['value']
        zone_start, zone_end = self.get_friction_zone(zone_index)
        for node_index in range(zone_start, zone_end + 1):
            self.set_friction_minor({'index': node_index, 'value': zone_value})
        self.logger.debug('Zone Ks minor set.')

    def set_zone_friction_major(self, ks):
        """Change major friction coefficient at zone.

        Uses :meth:`get_friction_zone` and :meth:`set_friction_major`.
        Does not use get_zone_friction_major but loop on set_friction_major
        as perturbation to ks is not additive.

        :param dict ks: Friction coefficient at zone ``{'index','value'}``
        """
        zone_index, zone_value = ks['index'], ks['value']
        zone_start, zone_end = self.get_friction_zone(zone_index)
        for node_index in range(zone_start, zone_end + 1):
            self.set_friction_major({'index': node_index, 'value': zone_value})
        self.logger.debug('Zone Ks major set.')

    def get_friction_minor(self, index=None):
        """
        Get minor friction coefficient at given node.

        Uses :meth:`Mascaret.get`.

        :param int index: node index
        :return: Minor friction coefficient
        :rtype: [list(]float[)]
        """
        k_s = [self.masc.get('Model.FricCoefMainCh', i)
               for i in range(self.model_size)]
        if index is not None:
            k_s = k_s[index]
        self.logger.debug('Ks minor value= {}.'.format(k_s))
        return k_s

    def get_friction_major(self, index=None):
        """Get major friction coefficient at given node.

        Uses :meth:`Mascaret.get`.

        :param int index: node index
        :return: Major friction coefficient
        :rtype: [list(]float[)]
        """
        k_s = [self.masc.get('Model.FricCoefFP', i)
               for i in range(self.model_size)]
        if index is not None:
            k_s = k_s[index]
        self.logger.debug('Ks major value= {}.'.format(k_s))
        return k_s

    def set_friction_minor(self, k_s):
        """Change minor friction coefficient.

        Uses :meth:`Mascaret.set`.
        Does not use get_friction minor because modif ks
         is not additive to old value

        :param dict k_s: Minor friction coefficient ``{'index','value'}``
        """
        self.logger.debug('Ks minor new value= {}'.format(k_s['value']))
        self.masc.set('Model.FricCoefMainCh', k_s['value'], k_s['index'])
        self.logger.debug('Ks minor changed.')

    def set_friction_major(self, k_s):
        """Change major friction coefficient.

        Uses :meth:`Mascaret.set`.
        Does not use get_friction major because modif
         ks is not additive to old value
        :param dict k_s: Major friction coefficient ``{'index','value'}``
        """
        self.logger.debug('Ks major new value= {}'.format(k_s['value']))
        self.masc.set('Model.FricCoefFP', k_s['value'], k_s['index'])
        self.logger.debug('Ks major changed.')

    def set_time_step(self, t_s):
        """Change time step.

        Uses: meth:`Mascaret.set`
        :param float t_s: time step
        """

        self.masc.set('Model.DT', t_s)
        self.time_step = t_s

    def perturbate_cross_sections(self, perturbation):
        """Change cross sections z by:
        * an additive term +dz everywhere
        * or by a Gaussian process instance

        Uses :meth:`Mascaret.get_var_size`, :meth:`Mascaret.get`
        and :meth:`Mascaret.set`.

        :param dict perturbation: Cross section perturbation
        """
        nb_secs, _, _ = self.masc.get_var_size('Model.CrossSection.RelAbs')
        nb_pts = []
        for sec in range(nb_secs):
            _, idx_nb_pts, _ = self.masc.get_var_size('Model.CrossSection.X',
                                                      sec)
            nb_pts.append(idx_nb_pts)
        print('nb of section', nb_secs)
        print('nb of pts per section', nb_pts)

        print('in Pertubate cross section', perturbation)
        shift_dz = np.zeros(nb_secs, float)
        if 'variation' in perturbation:
            shift_dz = shift_dz + perturbation['variation']
        print('shift_dz=', shift_dz)

        if 'section' in perturbation:
            sec = perturbation['section']
            self.logger.debug(
                'In Setter Cross Section, profil idx = {}'.format(sec))
            old_z_bot = self.cross_sections['bottom']['level'][sec]
            new_z_bot = old_z_bot + shift_dz[sec]
            self.masc.set('Model.CrossSection.Zbot', new_z_bot, sec)
            self.cross_sections['bottom']['level'][sec] = new_z_bot
            for pt in range(nb_pts[sec]):
                old_z = self.cross_sections['all']['level'][sec][pt]
                new_z = old_z + shift_dz[sec]
                self.masc.set('Model.CrossSection.Y', new_z, sec, pt)
                self.cross_sections['all']['level'][sec][pt] = new_z
                self.logger.info("Modified Cross Section: {} {}"
                                 .format(sec, self.cross_sections['all']
                                                                 ['level'][sec]
                                                                 [pt]))
        else:
            # Loop on number of section
            for sec in range(nb_secs):
                old_z_bot = self.cross_sections['bottom']['level'][sec]
                new_z_bot = old_z_bot + shift_dz[sec]
                self.masc.set('Model.CrossSection.Zbot', new_z_bot, sec)
                self.cross_sections['bottom']['level'][sec] = new_z_bot
                # Loop on number of point by section
                for pt in range(nb_pts[sec]):
                    old_z = self.cross_sections['all']['level'][sec][pt]
                    new_z = old_z + shift_dz[sec]
                    self.masc.set('Model.CrossSection.Y', new_z, sec, pt)
                    self.cross_sections['all']['level'][sec][pt] = new_z
                    print('sec,pt,new_z', sec, pt, new_z)

    def plot_water(self, xlab='Curvilinear abscissa (km)',
                   ylab1='Water elevation (m)',
                   ylab2='Upstream discharge (m3/s)',
                   title='Water evelation along the \
                   open-channel at final time',
                   output='WaterElevation'):
        """Plot the water elevation at final time

        :param str xlab: label s
        :param str ylab1: label y1
        :param str ylab2: label y2
        :param str title: title
        :param str output: output file
        """

        abscissa_km = copy.deepcopy(self.curvilinear_abscissa)
        abscissa_km[:] = [sk / 1000. for sk in self.curvilinear_abscissa]
        bathymetry = self.global_bathymetry
        elevation = self.global_elevation
        discharge = self.global_discharge
        fig, ax1 = plt.subplots()
        ax1.plot(abscissa_km, bathymetry, color='black')
        ax1.plot(abscissa_km, elevation, color='blue')
        ax1.fill_between(abscissa_km, bathymetry, elevation,
                         facecolor='blue', alpha=0.5)
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab1, color='blue')
        ax1.tick_params('y', colors='blue')
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax2 = ax1.twinx()
        ax2.plot(abscissa_km, discharge, color='red')
        ax2.set_ylabel(ylab2, color='red')
        ax2.tick_params('y', colors='red')
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        if 'output' in self.settings:
            hydraulic_state = self.local_hydraulic_state
            s_km = hydraulic_state['s'] / 1000.
            z_b = hydraulic_state['zb']
            z_val = hydraulic_state['z']
            q_val = hydraulic_state['q']
            plt.axvline(x=s_km, color='k', linestyle='dashed', linewidth=0.5)
            ax1.axhline(y=z_b, color='k', linestyle='dashed', linewidth=0.5)
            ax1.axhline(y=z_val, color='b', linestyle='dashed', linewidth=0.5)
            ax2.axhline(y=q_val, color='r', linestyle='dashed', linewidth=0.5)
            plt.title(title + '\n s={}km, zb={}m, z={}m, q={}m3/s'
                      .format(s_km, np.round(z_b, 2), np.round(z_val, 2),
                              int(q_val)))
        else:
            plt.title(title)
        fig.tight_layout()
        fig.savefig(os.path.join(self.paths['figures'], output + '.pdf'),
                    transparent=True, bbox_inches='tight')
        plt.close('all')

    def save(self, out_name=None, out_idx=0):
        np.savetxt(os.path.join(self.paths['output'],
                                'abscissa_' + out_name + str(out_idx) +
                                '.txt'),
                   self.curvilinear_abscissa)
        np.savetxt(os.path.join(self.paths['output'],
                                'elevation_' + out_name + str(out_idx) +
                                '.txt'),
                   self.global_elevation)
        np.savetxt(os.path.join(self.paths['output'],
                                'discharge_' + out_name + str(out_idx) +
                                '.txt'),
                   self.global_discharge)
        np.savetxt(os.path.join(self.paths['output'],
                                'bathymetry_' + out_name + str(out_idx) +
                                '.txt'),
                   self.global_bathymetry)
        if 'output' in self.settings:
            hydraulic_state = self.local_hydraulic_state
            np.savetxt(os.path.join(self.paths['output'],
                                    'local_abscissa_' + out_name +
                                    str(out_idx) + '.txt'),
                       [hydraulic_state['s']])
            np.savetxt(os.path.join(self.paths['output'],
                                    'local_elevation_' + out_name +
                                    str(out_idx) + '.txt'),
                       [hydraulic_state['z']])
            np.savetxt(os.path.join(self.paths['output'],
                                    'local_discharge_' + out_name +
                                    str(out_idx) + '.txt'),
                       [hydraulic_state['q']])
            np.savetxt(os.path.join(self.paths['output'],
                                    'local_bathymetry_' + out_name +
                                    str(out_idx) + '.txt'),
                       [hydraulic_state['zb']])

    def plot_bathymetry(self, xlab='Curvilinear abscissa (km)',
                        ylab='Bathymetry elevation (m)',
                        title='Bathymetry elevation along the open-channel',
                        output='BathymetryElevation'):
        """Plot the bathymetry elevation.

        :param str xlab: label s
        :param srt ylab: lable y
        :param str title: title
        :param str output: output file
        """
        abscissa_km = copy.deepcopy(self.curvilinear_abscissa)
        abscissa_km[:] = [sk / 1000. for sk in self.curvilinear_abscissa]
        plt.plot(abscissa_km, self.global_bathymetry, color='b',
                 marker='.', linestyle='-', markersize=1)
        if 'output' in self.settings:
            s = self.local_hydraulic_state['s'] / 1000.
            plt.axvline(x=s, color='k', linestyle='dashed')
        plt.ylabel(ylab)
        plt.xlabel(xlab)
        plt.title(title)
        plt.savefig(os.path.join(self.paths['figures'], output + '.pdf'),
                    transparent=True, bbox_inches='tight')
        plt.close('all')

    def plot_misc(self, plotlist, coords=None, color=None, labels=None,
                  xlabel='', ylabel='', name='plot', axis2=None, color2='g',
                  marker2='*', ylabel2=''):
        """Plot miscelaneous fields

        :plotlist: list of fields
        :coords: list of coordinates
        :color: color for each field
        :labels: labels for each field
        :ylabel: label for y axis
        :name:  name of plot
        :axis2: field for second axis
        :color2: color for second axis
        :marker2: marker for second axis
        :ylabel2: label for second y axis
        """
        if coords is None:
            abscissa_km = copy.deepcopy(self.curvilinear_abscissa)
            abscissa_km[:] = [sk / 1000. for sk in self.curvilinear_abscissa]
        else:
            abscissa_km = copy.deepcopy(coords)
        #            abscissa_km[:] = [sk/1000. for sk in coords]

        if color is None:
            color = ['k', 'b', 'r', 'g', 'c', 'm', 'y']

        fig, ax1 = plt.subplots()
        for cnt, var in enumerate(plotlist):
            if labels is None:
                ax1.plot(abscissa_km, var, color[cnt], linewidth=2.)
            else:
                ax1.plot(abscissa_km, var, color[cnt], linewidth=2.,
                         label=labels[cnt])
        ax1.set_ylabel(ylabel)
        ax1.set_xlabel('Curvilinear abscissa (km)')
        if xlabel is None:
            ax1.set_xlabel('Curvilinear abscissa (km)')
        else:
            ax1.set_xlabel(xlabel)
        ax1.grid()
        if labels is not None:
            plt.legend()

        if axis2 is not None:
            ax2 = ax1.twinx()
            ax2.plot(abscissa_km, axis2, '{}{}'.format(color2, marker2))
            ax2.set_ylabel(ylabel2, color=color2)
            ax2.tick_params('y', colors=color2)

        plt.savefig(os.path.join(self.paths['figures'], '{}.pdf'.format(name)),
                    format='pdf', transparent=True, bbox_inches='tight')
        plt.close('all')

    def get_info_init_file(self):
        """ Get information fromp the header of the initialization file.
        """
        in_file = open(os.path.join(self.paths['data'],
                                    self.settings['files']['lig']), "rt")
        contents = in_file.read().split()
        in_file.close()
        if 'IMAX' in contents:
            imax = int(contents[contents.index("IMAX") + 2])
        else:
            imax = int(contents[contents.index("IMAX=") + 1])
        if 'NBBIEF' in contents:
            nbbief = int(contents[contents.index("NBBIEF") + 2])
        else:
            nbbief = int(contents[contents.index("NBBIEF=") + 1])
        if 'I1,I2' in contents:
            i_1 = int(contents[contents.index("I1,I2") + 2])
            i_2 = int(contents[contents.index("I1,I2") + 3])
        else:
            i_1 = int(contents[contents.index("NBBIEF=") + 1])
            i_2 = int(contents[contents.index("NBBIEF=") + 2])
        x_val = [float(contents[i])
                 for i in range(contents.index("X") + 1, contents.index("Z"))]
        z_val = [float(contents[i])
                 for i in range(contents.index("Z") + 1, contents.index("Q"))]
        q_val = [float(contents[i])
                 for i in range(contents.index("Q") + 1,
                                contents.index("FIN"))]
        return {'IMAX': imax,
                'NBBIEF': nbbief,
                'I1': i_1,
                'I2': i_2,
                'X': x_val,
                'Z': z_val,
                'Q': q_val}

    def get_info_restart(self):
        """ Get information for dumping a restart file.
        """
        imax = self.model_size
        nbbief = self.masc.get_var_size('Model.Connect.FirstNdNum')[0]
        i_1 = []
        i_2 = []
        for i in range(nbbief):
            i_1.append(self.masc.get('Model.Connect.FirstNdNum', i))
            i_2.append(self.masc.get('Model.Connect.LastNdNum', i))
        return {'IMAX': imax,
                'NBBIEF': nbbief,
                'I1': i_1,
                'I2': i_2}

    def store_xcoord(self):
        """
        Compute the local relative abscissa as used in lig restart files
        """
        # Size informations
        imax = self.model_size
        nbbf = self.masc.get_var_size('Model.Connect.FirstNdNum')[0]
        # Local abscissa
        oribf = [self.masc.get('Model.Connect.FirstNdNum', i) - 1
                 for i in range(nbbf)]
        endbf = [self.masc.get('Model.Connect.LastNdNum', i)
                 for i in range(nbbf)]
        # Assign the bief number to each section (piecewise constant list)
        ibief = [ib + 0 * i for ib in range(nbbf) for i in range(oribf[ib],
                                                                 endbf[ib])]

        self.xcoord = []
        for i in range(imax):
            self.xcoord.append(
                self.masc.get('Model.X', i) -
                self.masc.get('Model.X', oribf[ibief[i]]) +
                self.masc.get('Model.CrossSection.RelAbs',
                              self.masc.get('Model.IDT', oribf[ibief[i]]) - 1))

        del ibief
        del oribf
        del endbf

    def save_lig_restart(self, out_file='RestartLigneEau.lig', k_s=None):
        """
        Save a lig restart files

        :param str out_file: name (and path) of the output restart file
        :param str k_s:       OPTIONAL : Ks other than None activate
            CF1 and CF2 storage
        """
        # Size informations
        imax = self.model_size
        nbbf = self.masc.get_var_size('Model.Connect.FirstNdNum')[0]
        # Formatting parameters
        nbcol = 5
        nbent = int((nbbf - 1) / nbcol) + 1
        nbfmt = int(imax / nbcol) * nbcol
        nbmod = imax % nbcol
        nblin = int(nbfmt / nbcol)
        # Filename
        with open(out_file, 'w') as lig:
            # Header
            lig.write('RESULTATS CALCUL, DATE :  '
                      '{}\n'.format(datetime.now().strftime('%d/%m/%Y %H:%M')))
            lig.write('FICHIER RESULTAT MASCARET\n')
            lig.write('------------------------------------------------\
                -----------------------\n')
            lig.write(' IMAX =%6i ' % imax)
            lig.write('NBBIEF=%5i\n' % nbbf)
            lig.write((nbent * ' ENTETE NON RELUE\n'))
            # X
            np.savetxt(lig,
                       np.asarray(self.xcoord[0:nbfmt]).reshape(nblin, nbcol),
                       fmt=(nbcol * '%13.2f'), header='X', comments=' ')
            if nbmod != 0:
                np.savetxt(lig,
                           np.asarray(self.xcoord[nbfmt:]).reshape(1, nbmod),
                           fmt=(nbmod * '%13.2f'))
            #  Z
            np.savetxt(lig,
                       np.asarray(self.global_elevation[0:nbfmt]).reshape(
                            nblin, nbcol),
                       fmt=(nbcol * '%13.3f'), header='Z', comments=' ')
            if nbmod != 0:
                np.savetxt(lig,
                           np.asarray(self.global_elevation[nbfmt:]).
                           reshape(1, nbmod),
                           fmt=(nbmod * '%13.3f'))
            # Q
            np.savetxt(lig,
                       np.asarray(self.global_discharge[0:nbfmt]).
                       reshape(nblin, nbcol),
                       fmt=(nbcol * '%13.3f'), header='Q', comments=' ')
            if nbmod != 0:
                np.savetxt(lig,
                           np.asarray(self.global_discharge[nbfmt:]).
                           reshape(1, nbmod),
                           fmt=(nbmod * '%13.3f'))

            # Optional Friction coefficients
            if k_s is not None:
                # CF1
                k_s = self.get_friction_minor()
                np.savetxt(lig,
                           np.asarray(k_s[0:nbfmt]).reshape(nblin, nbcol),
                           fmt=(nbcol * '%13.3f'), header='CF1', comments=' ')
                if nbmod != 0:
                    np.savetxt(lig,
                               np.asarray(k_s[nbfmt:]).reshape(1, nbmod),
                               fmt=(nbmod * '%13.3f'))

                # CF2
                k_s = self.get_friction_major()
                np.savetxt(lig,
                           np.asarray(k_s[0:nbfmt]).reshape(nblin, nbcol),
                           fmt=(nbcol * '%13.3f'), header='CF2', comments=' ')
                if nbmod != 0:
                    np.savetxt(lig,
                               np.asarray(k_s[nbfmt:]).reshape(1, nbmod),
                               fmt=(nbmod * '%13.3f'))

            # Footer
            lig.write(' FIN')
        lig.close()
