#!/usr/bin/env python
"""
Main programm of running coupling
=================
"""
import datetime
import glob
import json
import math
import os
import re
import subprocess
import tempfile
import time
from shutil import copy2, rmtree, move, copystat
from config import CFGS, update_config
from data_manip.extraction.telemac_file import TelemacFile
from execution.study import Study
from execution.telemac_cas import TelemacCas
from pretel.chop_selafin import ChopSelafin
from telapy.tools.driven_utils import shell_cmd, mpirun_cmd


def replace_all(txt, dico):
    """
    Replace several items

    @param txt (str) orginal text
    @param dico (dict) dictionary of expressions to replace

    @return (str) changed text
        """
    for i in dico:
        txt = txt.replace(i, dico[i])
    return txt


def file_path(ligne):
    """
    To complete the path file
    @param ligne (str) line
    @return (str) path file
    """
    if len(ligne) > 2:
        complet = ''
        for path in ligne:
            complet = os.path.join(complet, os.path.expandvars(path))
        config_file = complet
    else:
        config_file = os.path.expandvars(ligne[1])
    return config_file


def checkfile(cas_file, txt=None):
    """
    Check if file existes
    @param cas_file (str) name file
    @param txt (str)  text to print
    """
    if not os.path.exists(cas_file):
        if txt is not None:
            raise IOError(txt)
        else:
            raise IOError('\nFILE MISSING: {}. STOP. \n\n'.format(cas_file))


def file_mv(ori, dest, verbose=True):
    """
    Move file
    @param ori (str) path source
    @param dest (str) path target
    @param verbose (bool) display text
    """
    if os.path.exists(ori):
        move(ori, dest)
    else:
        if verbose:
            print('\nFILE MISSING {} FOR MOVE. \n\n'.format(ori))


class LongCplDriver:
    """
    Class coupling management
    """

    def __init__(self, case_path=".", coupling_def=None,
                 models_configs=None,
                 keep_exec=False):
        """
        Constructor
        @param case_path (str) path of model
        @param coupling_def (dict) coupling parameter
        @param models_configs (dict) dictionnary of submodel parameter
        @param keep_exec (bool) preserve the EXEC directory after execution
        """
        if not os.path.isdir(case_path):
            raise Exception('ERROR: CPLROOT dir {} does not exist\n'
                            '\nSTOP'.format(case_path))

        self.case_path = os.path.realpath(case_path)

        if coupling_def is None:
            raise NameError('\nDICO NOT FOUND {}. '
                            'PROGRAM STOP. \n\n'.format('coupling_def'))
        self.couplingdef = coupling_def

        if models_configs is None:
            raise NameError('\nDICO NOT FOUND {} . '
                            'PROGRAM STOP. \n\n'.format('models_configs'))
        self.list_dico_mod = models_configs

        self.keep_exec = keep_exec

        self.cas_file = None
        self.cas_file_tmp = None
        self.systelcfg = None
        self.usetelcfg = None
        self.dico_file = None
        self.configrun = None
        self.dico_file_model1d = {}

        self.param = {}
        self.param['dt_couplage'] = self.couplingdef["Coupling"]["TimeStep"]
        self.param['nb_mod_1d'] = len(self.couplingdef["1D"])
        self.param['nom_mod_1d'] = list(self.couplingdef["1D"].keys())
        self.param['dt_1d'] =  \
            self.couplingdef["1D"][self.param['nom_mod_1d'][0]]["TimeStep"]
        self.param['nb_mod_2d'] = len(self.couplingdef["2D"])
        self.param['nom_mod_2d'] = list(self.couplingdef["2D"].keys())
        self.param['dt_2d'] =  \
            self.couplingdef["2D"][self.param['nom_mod_2d'][0]]["TimeStep"]
        self.param['freq_res'] = \
            self.couplingdef["2D"][self.param['nom_mod_2d'][0]]["OutputFreq"] \
            / self.param['dt_couplage']

    def lec_telfile(self, dico):
        """
        Import data from the 2D model description dict

        @param dico (dict) Contains model information
        @return  (dict) , (dict):
            parameters dictionary of "cas" file  and "config" file
        """

        # DBG print("Read in the T2D cas file {}".format(dico['cas']))
        checkfile(dico['cas'])
        self.cas_file = dico['cas']

        if 'config_file' in dico:
            if dico['config_file'].lower().find('${systelcfg}') > -1:
                dico_repl = {'${systelcfg}': str(os.environ['SYSTELCFG']),
                             '${SYSTELCFG}': str(os.environ['SYSTELCFG'])}
                self.systelcfg = replace_all(dico['config_file'], dico_repl)
            else:
                self.systelcfg = dico['config_file']
        else:
            self.systelcfg = str(os.environ['SYSTELCFG'])

        # DBG print("Using the {} config file".format(self.systelcfg))
        checkfile(self.systelcfg)

        if 'config_option' in dico:
            if dico['config_option'].lower().find('${usetelcfg}') > -1:
                dico_repl = {'${usetelcfg}': str(os.environ['USETELCFG']),
                             '${USETELCFG}': str(os.environ['USETELCFG'])}
                self.usetelcfg = replace_all(dico['config_option'], dico_repl)
            else:
                self.usetelcfg = dico['config_option']
        else:
            self.usetelcfg = str(os.environ['USETELCFG'])

        # DBG print("Using the {} config option".format(self.usetelcfg))

    def suffixe_fichier(self, liste):
        """
        Delete file suffix
        @param liste (list)
        """
        if len(liste) == 1:
            id_proc_2d = liste[0]
            self.param['nb_proc_2D'] = 1
        elif len(liste) == 2:
            id_proc_2d = liste[0]
            self.param['nb_proc_2D'] = liste[1]
        else:
            id_proc_2d = 1
            self.param['nb_proc_2D'] = 1

        id_proc = "{:05}".format(id_proc_2d)
        nb_proc = "{:05}".format(self.param['nb_proc_2D'])

        if len(liste) == 1:
            print(id_proc)
        elif len(liste) == 2:
            print(int(nb_proc) - int(id_proc))

    def exec_creation(self):
        """
        Creation EXEC represitory
        """
        if os.path.isdir(os.path.join(self.case_path, 'EXEC')):
            rmtree(os.path.join(self.case_path, 'EXEC'))
            os.mkdir(os.path.join(self.case_path, 'EXEC'))
        else:
            os.mkdir(os.path.join(self.case_path, 'EXEC'))

    def preparation_model(self):
        """
        Preparation of model
        1. copy file for coupling
        2. Modification parameters file
        """
        # # ~~~ Identification of TELAMAC files ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        tel_file = "config_{}".format(self.param['nom_mod_2d'][0])

        if tel_file in self.list_dico_mod.keys():
            dico = self.list_dico_mod[tel_file]
        else:
            raise NameError('\nDICO NOT FOUND {}. '
                            'PROGRAM STOP. \n\n'.format(tel_file))
        self.lec_telfile(dico['files'])

        tstart = self.param['startTime']

        # configuration define
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        root_dir = os.path.expandvars('$HOMETEL')
        python_dir = os.path.join(root_dir, 'scripts', 'python3')
        CFGS.parse_cfg_file(self.systelcfg, self.usetelcfg, root_dir,
                            python_dir)
        CFGS.compute_modules_info()
        CFGS.compute_system_info()
        CFGS.compute_partel_info()
        CFGS.compute_mpi_info()

        #
        # # ~~~ Modification of .cas files et config files in function
        # #      to 2D run type  ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.param['Run_type_2D'].lower() == 'sequentiel':
            self.param['nb_proc_2D'] = 1

        # # ~~~ Creating symbolic link for initial water line if reprise ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        folder_cas = os.path.dirname(self.cas_file)

        if self.param['reprise']:
            file_tmp = os.path.join(folder_cas, 'WaterLineInit_{}.slf'
                                    .format(tstart))
            if not os.path.exists(file_tmp):
                raise Exception("ERROR. INTIAL WATER LINE"
                                "OF 2D MODEL IS MISSED."
                                " STOP.")
            else:
                path_tmp2 = os.path.join(folder_cas, "WaterLineInit_in.slf")
                if os.path.exists(path_tmp2):
                    os.remove(path_tmp2)
                path_ori = 'WaterLineInit_{}.slf'.format(tstart)
                os.symlink(path_ori, path_tmp2)

        #
        i_b = 0
        while i_b < self.param['nb_mod_2d']:
            file_tmp = os.path.join(folder_cas, 'bc2D_restart_{}_{}.json'
                                    .format(self.param['nom_mod_2d'][i_b],
                                            tstart))
            if os.path.exists(file_tmp):
                copy2(file_tmp, '.')
            i_b += 1

        # print("~~~~~ CREATION TEMPORARY DATA FILE OF TELEMAC ~~~~~")
        os.chdir(os.path.realpath(os.path.dirname(self.cas_file)))
        my_study = Study(os.path.basename(self.cas_file), 'telemac2d',
                         self.cas_file_tmp)
        os.chdir(self.cas_file_tmp)
        self.dico_file = my_study.dico_file
        # creation  split files
        my_study.fill_working_dir()
        my_study.ncsize = self.param['nb_proc_2D']
        my_study.generate_mpi_files()

        #
        # ~~~ Test the consistence between the last time
        # of restart file and specific time ~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if self.param['reprise']:
            res = TelemacFile('T2DPRE')
            last_time_t2dpre = float(res.times[-1])
            res.close()
            if last_time_t2dpre != tstart:
                raise Exception(
                    "\nERROR, THE RESTART TIME OF PARAMETERS FILES {} "
                    "AND OF RESTART {} AREN'T THE SAME."
                    "PROGRAM STOP.\n\n".format(tstart, last_time_t2dpre))

        # # ~~~ Check data of Mascaret  ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        i_b = 0
        while i_b < self.param['nb_mod_1d']:
            name_dict = "config_{}".format(self.param['nom_mod_1d'][i_b])
            if name_dict in self.list_dico_mod.keys():
                dico = self.list_dico_mod[name_dict]
                self.dico_file_model1d[self.param['nom_mod_1d'][i_b]] = \
                    dico['files']
            i_b += 1
        #
        #
        # # ~~~ Modification .cas file from coupling data ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #
        # # ~~~ 2D.cas  file~~~~~
        # print('~~~~~~~~~~Modification of .cas file~~~~~~~~~~\n')
        periode_sorties = \
            int(self.param['freq_res'] * self.param['dt_couplage']
                / self.param['dt_2d'])
        filecas = os.path.join(self.cas_file_tmp, 'T2DCAS')
        cas = TelemacCas(filecas, self.dico_file, check_files=False)

        cas.set('TIME STEP', self.param['dt_2d'])
        cas.set('GRAPHIC PRINTOUT PERIOD', periode_sorties)
        cas.set('LISTING PRINTOUT PERIOD', periode_sorties)
        cas.set('GEOMETRY FILE', 'T2DGEO')
        cas.set('BOUNDARY CONDITIONS FILE', 'T2DCLI')
        cas.set('RESULTS FILE', 'T2DRES')
        if 'LIQUID BOUNDARIES FILE' in cas.values.keys():
            cas.set('LIQUID BOUNDARIES FILE', 'T2DIMP')

        if 'PREVIOUS COMPUTATION FILE' in cas.values.keys():
            cas.set('PREVIOUS COMPUTATION FILE', 'T2DPRE')

        if self.param['reprise']:
            cas.set('COMPUTATION CONTINUED', True)
        else:
            cas.set('COMPUTATION CONTINUED', False)

        cas.write(filecas)

    def __call__(self, config_run=None):
        """
        Function controlling the coupling
        """
        # ---------------------------------------------------------
        # ~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~ CONFIGURATION ~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~

        if config_run is None:
            raise NameError('\nDICO NOT FOUND {}. '
                            'PROGRAM STOP. \n\n'.format('config_run'))
        self.configrun = config_run

        dicorun = self.configrun["Run"]
        ref_date = datetime.datetime.strptime(dicorun["RefDate"],
                                              '%d/%m/%Y %H:%M:%S')
        start_date = datetime.datetime.strptime(dicorun["StartDate"],
                                                '%d/%m/%Y %H:%M:%S')
        end_date = datetime.datetime.strptime(dicorun["EndDate"],
                                              '%d/%m/%Y %H:%M:%S')
        if "SingleExecDuration" in dicorun:

            iso8601_duration_re = re.compile(
                r'^(?P<sign>[-+]?)'
                r'(?:(?P<days>\d+(.\d+)?)d)?'
                r'(?:(?P<hours>\d+(.\d+)?)h)?'
                r'(?:(?P<minutes>\d+(.\d+)?)m)?'
                r'(?:(?P<seconds>\d+(.\d+)?)s)?'
                r'$'
            )

            postgres_interval_re = re.compile(
                r'^'
                r'(?:(?P<days>-?\d+) (days? ?))?'
                r'(?:(?P<sign>[-+])?'
                r'(?P<hours>\d+):'
                r'(?P<minutes>\d\d):'
                r'(?P<seconds>\d\d)'
                r')?$'
            )

            match = (postgres_interval_re.match(dicorun["SingleExecDuration"])
                     or
                     iso8601_duration_re.match(dicorun["SingleExecDuration"]))
            k_w = match.groupdict()
            days = datetime.timedelta(float(k_w.pop('days', 0) or 0))
            sign = -1 if k_w.pop('sign', '+') == '-' else 1
            k_w = {k: float(v) for k, v in k_w.items() if v is not None}
            length = days + sign * datetime.timedelta(**k_w)
        else:
            length = end_date - start_date

        self.param['startTime'] = int((start_date - ref_date).total_seconds())
        self.param['endTime'] = int((end_date - ref_date).total_seconds())
        self.param['singleRun'] = int(length.total_seconds())
        if "SaveCheckPoints" in dicorun:
            self.param['chkPts'] = dicorun["SaveCheckPoints"].lower() == "yes"
        else:
            self.param['chkPts'] = False

        self.param['reprise'] = dicorun["RestartFromFile"].lower() == "yes"

        if self.configrun["2D"][self.param['nom_mod_2d'][0]]["Parallel"] \
                .lower() == 'no':
            self.param['Run_type_2D'] = 'sequentiel'
            self.param['nb_proc_2D'] = 1
        else:
            self.param['Run_type_2D'] = 'parallele'
            self.param['nb_proc_2D'] = \
                self.configrun["2D"][self.param['nom_mod_2d'][0]]["NbProc"]

        nbrun = \
            math.ceil(float(self.param['endTime'] - self.param['startTime']) /
                      float(self.param['singleRun']))
        print("\n+--------------------------------------------------------+"
              "\n|------- LONGITUDINAL MASCARET TELEMAC2D COUPLING -------|"
              "\n+--------------------------------------------------------+"
              "\n|")
        print('|  Coupling Method:', self.couplingdef["Coupling"]["Method"])
        print('|  Max It. nr.:    ', self.couplingdef["Coupling"]["MaxIter"])
        print('|  Coupling TStep: ', self.couplingdef["Coupling"]["TimeStep"])
        print('|')
        print('|  Initial time:   ', dicorun["StartDate"])
        print('|  End time:       ', dicorun["EndDate"])
        if "SingleExecDuration" in dicorun:
            print('|  Split Run every:',
                  dicorun["SingleExecDuration"], "({} exec[s])"
                  .format(nbrun))
        else:
            print('|  Run not splitted ({} exec)'.format(nbrun))
        print('|  Restarted run:  ', dicorun["RestartFromFile"])
        if self.param['chkPts']:
            print('|  Checkpointed:    yes')
        print('|')
        print('|  1D models:      ', self.param['nom_mod_1d'])
        print('|  2D models:      ', self.param['nom_mod_2d'])
        print('|')
        print('|  Interfaces:      ')
        for itf in self.couplingdef["Interfaces"]:
            print('|    {:<12.12s} ({:<10.10s}) {:<10.10s} of {:<12.12s}'
                  .format("{}:{}".format(itf["Id1D"], itf["IdExtr1D"]),
                          itf["Condition1D"], itf["1DPosition"],
                          "{}:{}".format(itf["Id2D"], itf["LiqBdry2D"])))
        print("+--------------------------------------------------------+\n")

        # ~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~ PRELIMINAIRES ~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~

        start = time.time()
        previous_path = os.getcwd()
        os.chdir(self.case_path)

        self.exec_creation()
        os.chdir(os.path.join(self.case_path, 'EXEC'))
        self.cas_file_tmp = os.getcwd()
        self.preparation_model()
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # # ~~~ LANCEMENT DU MODELE COUPLE ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              "~~~~~   RUN THE COUPLED MODEL    ~~~~~\n"
              "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        tdeb = self.param['startTime']
        freq = self.param['singleRun']

        tfinal = tdeb + math.ceil((self.param['endTime'] - tdeb) / freq) * freq

        while tdeb < tfinal:
            tstart = tdeb
            tstop = tdeb + freq

            #
            #    # ~~~ MAJ des fichiers .cas ~~~~~
            #    # ~~~~~~
            #    # ~~~ Fichier .cas 2D ~~~~~

            filename = os.path.join(self.cas_file_tmp, 'T2DCAS')
            cas = TelemacCas(filename, self.dico_file, check_files=False)
            cas.set('DURATION', float(freq))

            if tdeb > self.param['startTime']:
                self.param['reprise'] = True
                cas.set('SUITE DE CALCUL', True)
            cas.write(filename)

            current_run = {"Run": {}}
            current_run["Run"]["StartTime"] = tstart
            current_run["Run"]["RunLength"] = freq
            current_run["Run"]["CplTSteps"] = \
                math.ceil(self.param['singleRun'] / self.param['dt_couplage'])
            current_run["Run"]["Restarted"] = self.param['reprise']
            current_run["Run"]["InitialRun"] = bool(tdeb == 0.)
            current_run["1D"] = {}
            i_b = 0
            while i_b < self.param['nb_mod_1d']:
                name = self.param['nom_mod_1d'][i_b]
                current_run["1D"][name] = {}
                folder_lig = \
                    os.path.dirname(self.dico_file_model1d[name]["lig"])
                file = "WaterLine_{}_{}.lig"\
                    .format(name, ('%.6f' % tstart).rstrip('0').rstrip('.'))
                file_lig = os.path.join(folder_lig, file)
                if not os.path.exists(file_lig):
                    raise IOError(
                        "ERROR: MISSING WATERLINE FILE {}".format(file_lig))
                current_run["1D"][name]["WaterLineFile"] = file_lig

                file = "bc1D_restart_{}_{}.json"\
                    .format(name, ('%.6f' % tstart).rstrip('0').rstrip('.'))
                file_bc = os.path.join(folder_lig, file)
                if os.path.exists(file_bc):
                    copy2(file_bc, os.path.join(self.case_path, 'EXEC'))

                current_run['config_{}'
                            .format(name)] =\
                    self.list_dico_mod['config_{}'.format(name)]

                i_b += 1
            current_run['coupling_def'] = self.couplingdef
            with open('CurrentRunDef.json', 'w') as f_p:
                json.dump(current_run, f_p, indent=4)
            #
            #    # ~~~ Creation of results directory from simulate time ~~~~~
            try:
                dirpath = 'COUPLING_FROM_{}'.format(tdeb)
                dirpath = os.path.join(self.case_path, 'Results', dirpath)
                os.makedirs(dirpath)
            except OSError:
                pass
            i_b = 0
            while i_b < self.param['nb_mod_1d']:
                pattern = 'listing .*$'
                filedst = os.path.join(self.case_path, "Results",
                                       "COUPLING_FROM_{}".format(tdeb),
                                       "ResultatsListing_{}.lis"
                                       .format(self.param['nom_mod_1d'][i_b]))
                repl = 'listing ' + filedst

                key = self.param['nom_mod_1d'][i_b]
                pattern_compiled = re.compile(pattern)
                line = self.dico_file_model1d[key]["listing"]
                self.dico_file_model1d[key]["listing"] \
                    = pattern_compiled.sub(repl, line)

                i_b += 1
            #    # ~~~  Suppression of old listing files of MASCARET~~~~~
            #    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            globpattern = os.path.join(self.case_path,
                                       'Results',
                                       'COUPLING_FROM_{}'.format(tdeb),
                                       'Res')
            globpattern = globpattern + "*.lis"
            for file in glob.glob(globpattern):
                os.remove(file)

            #    # ~~~ COMPUTE THE COUPLING ~~~~~
            #    # ~~~~~~~~~~~~~~~~~~~~~~~
            listing_file = os.path.basename(self.cas_file)
            listing_file = '{}.stdout'.format(listing_file)
            print("Running the coupled model between init time {} "
                  "and end time {}".format(tstart, tstop))

            cmd = mpirun_cmd()
            cmd = cmd.replace('<ncsize>', str(int(self.param['nb_mod_1d'] +
                              self.param['nb_proc_2D'])))

            launch_py = os.path.join(os.environ['HOMETEL'], 'scripts',
                                     'python3', 'run_cpl.py')
            launch_exe = "{0} launcher --n1d {1} > {2}"\
                .format(launch_py, self.param['nb_mod_1d'], listing_file)

            cmd = cmd.replace('<exename>', launch_exe)

            _, return_code = shell_cmd(cmd)

            if return_code != 0:
                raise Exception("\nERROR IN COUPLING RUN.\n"
                                "THE COMMAND IS : {} \n"
                                " PROGRAM STOP.\n".format(cmd))

            #
            #    # ~~~ Suppression of  MASCARET listing files~~~~~
            #    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            globpattern = os.path.join(self.case_path, 'Results',
                                       'COUPLING_FROM_{}'
                                       .format(tdeb), 'Res')
            globpattern = globpattern + "*.lis"
            for file in glob.glob(globpattern):
                os.remove(file)
            #
            #    # ~~~ Creation reprise files ~~~~~
            #    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #
            #    # ~~~ Creation 1D reprise files ~~~~~
            i_b = 0
            while i_b < self.param['nb_mod_1d']:
                folder_lig = os.path.dirname(
                    self.dico_file_model1d[
                        self.param['nom_mod_1d'][i_b]]["lig"])
                name_file = 'WaterLine_{}_*.lig'.format(
                    self.param['nom_mod_1d'][i_b])
                for file in glob.glob(name_file):
                    file_mv(file, os.path.join(folder_lig, file))
                if self.param['chkPts'] or tdeb + freq == tfinal:
                    name_file = 'bc1D_restart_{}_{}.json' \
                        .format(self.param['nom_mod_1d'][i_b], tdeb + freq)
                    copy2(name_file, os.path.join(folder_lig, name_file))
                i_b += 1
            if tdeb + freq < tfinal:
                #    # ~~~ Creation 2D reprise files ~~~~~
                chop_step = int(freq / (self.param['freq_res'] *
                                self.param['dt_couplage']) + 1)

                if self.param['Run_type_2D'].lower() == 'parallele':
                    prepath_para = os.path.join(self.cas_file_tmp, 'T2DPRE')
                    globpattern = prepath_para + "*"
                    for file in glob.glob(globpattern):
                        os.remove(file)

                slf = ChopSelafin('T2DRES', (chop_step, chop_step, chop_step))
                slf.put_content('T2DPRE')

            #    # ~~~ Creation result files ~~~~~
            #    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #
            #    # ~~~  .opt files~~~~~
            i_b = 0
            while i_b < self.param['nb_mod_1d']:
                filesrc = "study_par_{}".format(self.param['nom_mod_1d'][i_b])
                filesrc = os.path.join(filesrc, "output",
                                       "ResultatsOpthyca_{0}.opt"
                                       .format(self.param['nom_mod_1d'][i_b]))
                filedst = os.path.join(self.case_path,
                                       "Results", "COUPLING_FROM_{}"
                                       .format(tdeb),
                                       "ResultatsOpthyca_{}.opt"
                                       .format(self.param['nom_mod_1d'][i_b]))
                file_mv(filesrc, filedst)

                filesrc = "study_par_{}".format(self.param['nom_mod_1d'][i_b])
                filesrc = \
                    os.path.join(filesrc, "output",
                                 self.dico_file_model1d[
                                    self.param['nom_mod_1d'][i_b]]['listing'])
                filedst = \
                    os.path.join(self.case_path, "Results",
                                 "COUPLING_FROM_{}".format(tdeb),
                                 self.dico_file_model1d[
                                    self.param['nom_mod_1d'][i_b]]['listing'])
                file_mv(filesrc, filedst, verbose=False)
                i_b += 1

            filedst = os.path.join(self.case_path, "Results",
                                   "COUPLING_FROM_{}"
                                   .format(tdeb), listing_file)
            file_mv(listing_file, filedst)

            #    # ~~~  T2D .slf FILES  ~~~~~
            result_file = os.path.basename(self.cas_file)
            result_file = '{}_Resultats.slf'.format(result_file)
            restart_dir = os.path.dirname(self.cas_file)
            restart_file = os.path.join(restart_dir,
                                        'WaterLineInit_{}.slf'
                                        .format(tdeb + freq))

            filedst = os.path.join(self.case_path,
                                   "Results", "COUPLING_FROM_{}".format(tdeb),
                                   result_file)
            if self.param['chkPts'] or tdeb + freq == tfinal:
                copy2("T2DRES", restart_file)
            file_mv("T2DRES", filedst)

            if self.couplingdef["Coupling"]["Method"].lower() == \
                "additiveschwarz" and \
                    (self.param['chkPts'] or tdeb + freq == tfinal):
                restart_dir = os.path.dirname(self.cas_file)
                i_b = 0
                while i_b < self.param['nb_mod_2d']:
                    restart_file = 'bc2D_restart_{}_{}.json'. \
                        format(self.param['nom_mod_2d'][i_b], tdeb + freq)
                    copy2(restart_file,
                          os.path.join(restart_dir, restart_file))
                    i_b += 1

            filedst = os.path.join(
                self.case_path, "Results", "COUPLING_FROM_{}".format(tdeb),
                "Convergence_criteria.out")
            file_mv("Convergence_criteria.out", filedst)
            filedst = os.path.join(
                self.case_path, "Results", "COUPLING_FROM_{}".format(tdeb),
                "Convergence_criteria.csv")
            file_mv("Convergence_criteria.csv", filedst)
            tdeb = tdeb + freq

            if os.path.isdir(os.path.join(self.case_path, 'EXEC',
                             '__pycache__')):
                rmtree(os.path.join(self.case_path, 'EXEC', '__pycache__'))

        # # ~~~~~~~~~~~~~~~~~~~~~~
        # # ~~~ FINALISATION ~~~~~
        # # ~~~~~~~~~~~~~~~~~~~~~~
        #
        # # ~~~ RETURN INITIAL PATH ~~~~~
        os.chdir("..")
        # cd ..
        end = time.time()
        print("My work is done. Coupled job lasted : {} s \n"
              .format(end - start))
        os.chdir(previous_path)

    def run(self):
        """
        Alternative call method
        """
        self()

    def __del__(self):
        """Delete an instance of coupling driver"""
        if not self.keep_exec:
            rmtree(os.path.join(self.case_path, 'EXEC'))
