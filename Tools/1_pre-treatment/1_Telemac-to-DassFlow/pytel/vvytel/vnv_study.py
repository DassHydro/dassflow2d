"""
File containing the VnvStudy class that is used to run validation
"""
from os import path, chdir, makedirs, remove, environ
from collections import OrderedDict
from abc import ABC, abstractmethod
import subprocess
import time
from config import CFGS
from execution.run_cas import run_local_cas, run_hpc_cas
from execution.study import Study
from execution.get import get_hpc_cmd
from utils.exceptions import TelemacException
from utils.files import is_newer, put_file_content
from vvytel.vnv_tools import compute_norm, compute_diff, check_compatibility
from data_manip.extraction.telemac_file import TelemacFile
from data_manip.formats.mascaret_file import MascaretFile

class AbstractVnvStudy(ABC):
    """
    Defines all the actions for validation
    """
    def __init__(self, name, file_name, options):
        """
        Init function

        @param name (str) Name of the study
        @param file_name (str) Path to the validation file
        @param options (ArgumentParser) script options (from
        validate_telemac.py)
        """
        # Tag
        self.tag = []
        # Rank definition
        self.rank = 0
        # Files that will be created by the pre step (to be removed when doing a --clean)
        self.temporary_files = []

        # Name of the vnvStudy object
        self.name = name

        self.case_dir = path.dirname(file_name)
        self.file_name = file_name

        self.studies = OrderedDict()
        self.commands = OrderedDict()
        self.options = options
        self.action_time = OrderedDict()
        # Walltime for cluster run (default the one from options)
        self._default_walltime = options.walltime
        self.walltime = options.walltime
        self.listing = False

        self._init()

    def get_vnv_working_dir(self, name):
        """
        Returns the path of vnv working directory for a given study name

        @param name (str) Name of the vnv action

        @returns (str) The path to the vnv working directory
        """
        cfg = CFGS.configs[CFGS.cfgname]
        # The working dir is in
        # valdation folder/case_dir/val_name/vnv_name/<USETELCFG>
        vnv_working_dir = path.join(cfg['val_root'],
                                    self.case_dir,
                                    self.name,
                                    name,
                                    CFGS.cfgname)

        return vnv_working_dir

    def build_vnv_working_dir(self, name):
        """
        Builds the vnv working directory for a given action name

        @param name (str) Name of the vnv action

        @returns (str) The path to the vnv working directory
        """
        vnv_working_dir = self.get_vnv_working_dir(name)

        if not path.exists(vnv_working_dir):
            makedirs(vnv_working_dir)

        return vnv_working_dir

    def add_study(self, name, module, steering_file, cas=None):
        """
        Adding a study to the vnv class

        @param name (str) Name of the vnv study
        @param module (str) Name of the telemac-mascaret module
        @param steering_file (str) Name of the steering file
        @param cas (TelemacCas) If present create steering from data in cas and
        writing it as steering_file
        """
        if name in zip(self.studies.keys(), self.commands.keys()):
            raise TelemacException(\
                    "Study name ({}) already used pick another one\n"
                    "Names used: {}".format(name, self.studies.keys()))

        # Building vnv working dir
        vnv_working_dir = self.build_vnv_working_dir(name)

        # Path of the steering file in the working dir
        vnv_steering_file = path.join(\
                vnv_working_dir,
                path.basename(steering_file))

        # Special case were we give a TelemacCas structure instead of the
        # steering file
        if cas is not None:
            # Temporary creates a steering file prefixed with configuration for
            # when multiple configs are run at the same time
            tmp_steering_file = CFGS.cfgname+'_'+steering_file
            # Creating temporary steering file (to do the copy of the files)
            cas.write(tmp_steering_file)
            # Writing the one in the vnv working dir
            cas.write(vnv_steering_file)
            copy_cas = False
        else:
            tmp_steering_file = steering_file
            copy_cas = True

        # Temporary study just to copy the files
        local_study = Study(tmp_steering_file, module, '')

        local_study.copy_files(vnv_working_dir,
                               verbose=self.options.verbose,
                               copy_cas_file=copy_cas)

        del local_study

        working_dir = ''
        if module != 'mascaret':
            working_dir = path.join(vnv_working_dir, 'tmp_'+name)
        self.studies[name] = Study(vnv_steering_file, module, working_dir)
        if cas is not None:
            remove(tmp_steering_file)

    def add_command(self, name, command, hpc=False):
        """
        Adding a bash command to run
        """
        if name in zip(self.commands.keys(), self.studies.keys()):
            raise TelemacException(\
                    "Command name ({}) already used pick another one\n"
                    "Names used: {}".format(name, self.commands.keys()))

        self.commands[name] = {'cmd':command, 'hpc':hpc}

    def _run_case(self, name, study):
        """
        Running a case

        @param name (str) name of the vnv study
        @param study (study) study to run
        """
        #TODO: Add a silent option to the run of a case (at least for the
        # listing of telemac)
        vnv_working_dir = self.build_vnv_working_dir(name)

        chdir(vnv_working_dir)
        # Checking if the case is up to date then doing nothing
        # Criteria: one of the in_files is newer that the out_file
        # Path to the orignal case file (the one in the example folder)
        # This is to handles cases where that file does not exists because the
        ori_cas_file = path.join(self.case_dir,
                                 path.basename(study.steering_file))
        # steering file was generated
        cas_file = path.join(vnv_working_dir,
                             path.basename(study.steering_file))
        # For mascaret as it does not have a steering file in Telemac format we
        # always run it
        if study.code_name == 'mascaret':
            run = True
        else:
            # Take the first output file that is not an input file as well
            out_key = ''
            for key in study.cas.out_files:
                # Skipping PARAL files as the check for existence will fail
                # (the files names are change to add processor number)
                if 'PARAL' in study.cas.out_files[key]:
                    continue
                if out_key not in study.cas.in_files:
                    out_key = key
                    break

            if out_key == '':
                raise TelemacException(\
                        'Could not found a output file in your steering file'+
                        '\n{}'.format(study.cas.file_name))

            out_file = path.join(vnv_working_dir, study.cas.get(out_key))
            run = not path.exists(out_file)
            # checks if the steering file is newer than a result
            if path.exists(ori_cas_file) and is_newer(out_file, cas_file):
                run = True
            else:
                # Checks if any input files is newer than the output file
                for key in study.cas.in_files:
                    ffile = study.cas.get(key)
                    vnv_file = path.join(vnv_working_dir, ffile)
                    if is_newer(out_file, vnv_file):
                        run = True
                        break

        if run:
            print("  ~> {}: running on {} cores".format(name, study.ncsize))
            if self.listing:
                self.options.sortie_file = True
            # Forcing ncsize
            self.options.ncsize = study.ncsize
            # study specifiv walltime
            if isinstance(self.walltime, dict):
                if name in self.walltime:
                    self.options.walltime = self.walltime[name]
                else:
                    self.options.walltime = self._default_walltime
            else:
                self.options.walltime = self.walltime
            self.options.jobname = "{}__{}".format(path.basename(self.name), name)
            if study.cfg['HPC'] == {} or self.options.mpi:
                run_local_cas(study, self.options)
            else:
                run_hpc_cas(study, self.options)
        else:
            print("  ~> {}: Nothing to do (up-to-date)".format(name))

    def get_study_file(self, string):
        """
        Returns the file associated to a string containing 'name:submit' where
        name is the name in self.studies and submit is the name of the file in
        the telemac-mascaret temporary folder

        @param string (str) Name to convert

        @returns (str) Path to the file
        """
        name, submit = string.split(':')
        if name not in self.studies:
            raise TelemacException(\
                    "Unknown study name: {}\nKnown ones: {}"\
                    .format(name, self.studies.keys()))

        study = self.studies[name]

        file_name = None
        try:
            file_name = study.cas.get_file_from_submit(submit)
        except TelemacException:
            # Not in main steering file try cpl ones
            for cas in study.cpl_cases.values():
                try:
                    file_name = cas.get_file_from_submit(submit)
                except TelemacException:
                    continue
                # If we reach here the file was found
                break
            if file_name is None:
                raise TelemacException(\
                        "{} was not found in any of the steering files"\
                        .format(submit))
        # In case we have a file for each processor (PARAL in submit in
        # dictionary) get_file_from_submit will return a list
        if isinstance(file_name, list):
            vnv_wkg_dir = self.get_vnv_working_dir(name)
            return [path.join(vnv_wkg_dir, ffile) for ffile in file_name]

        return path.join(self.get_vnv_working_dir(name), file_name)

    def get_study_res(self, name=None, module="T2D",
                      whitelist=None, load_bnd=False):
        """
        Returns the telemac file result associated to a string containing
        'name:submit' where name is the name in self.studies and submit is
        the name of the file in the telemac-mascaret temporary folder. If
        name is None, returns a list of all results.

        @param name (str) Name of the study file (default: None)
        @param module (str) module of the res file (default: 'T2D')
        @param whitelist (list) select only result files that contains this
        element. Example of white list to select only sequential results:
        ['seq']. Strings composing the name must be separated by '_'. Return
        only the first part of the slitted name as result label (default: None)
        @param load_bnd (bool) load boundaries in telemac result file

        @returns (str) Path to the file
        """
        # Specific treatment of Mascaret file
        if module == 'mascaret':
            name1, submit = name.split(':')
            masc_file = path.join(self.get_vnv_working_dir(name1), submit)
            res = MascaretFile(masc_file)
            return res, name.upper()

        # Extract only one result file
        if name is not None:
            if ':' in name:
                name, submit = name.split(':')
            else:
                submit = module + 'RES'

            if name not in self.studies:
                raise TelemacException(\
                        "Unknown study name: {}\nKnown ones: {}"\
                        .format(name, self.studies.keys()))

            tel_file = self.get_study_file(name + ':' + submit)

            if load_bnd:
                bnd_file = self.get_study_file(name + ':' + \
                                               module + 'CLI')
                res = TelemacFile(tel_file, bnd_file=bnd_file)
            else:
                res = TelemacFile(tel_file)

            return res, name.upper()

        # Extract multiple result files and return a list
        res_list = []
        res_labels = []
        for nname in self.studies:
            # Apply white list
            if whitelist is not None:
                for white_str in whitelist:
                    if white_str in nname:
                        tel_file = self.get_study_file(nname + ':' + \
                                                       module + 'RES')
                        if load_bnd:
                            bnd_file = self.get_study_file(nname + ':' + \
                                                           module + 'CLI')
                            res_list.append(TelemacFile(tel_file,
                                                        bnd_file=bnd_file))
                        else:
                            res_list.append(TelemacFile(tel_file))

                        res_labels.append(nname.split('_')[0].upper())
            else:
                tel_file = self.get_study_file(nname + ':' + \
                                               module + 'RES')
                if load_bnd:
                    bnd_file = self.get_study_file(nname + ':' + \
                                                   module + 'CLI')
                    res_list.append(TelemacFile(tel_file,
                                                bnd_file=bnd_file))
                else:
                    res_list.append(TelemacFile(tel_file))

                res_labels.append(nname.upper())

        return res_list, res_labels

    def compute_errors(self,
                       res1=None, res2=None,
                       var1='all', var2='all', record=-1,
                       var3=None, var4=None,
                       relative=False, norm='linf', mass=None,
                       check_name=True):
        """
        Compute the norm of the difference between each variable of file1
        and file2 or between variables var1 and var2. Both result files must
        have the same mesh, the same variables and the same records.
        var1 and var2 can either be:
        - Name of a variable
        - 'all', all variables in file1 and file2 are checked

        @param res1 (TelemacFile) first telemac-mascaret file
        @param res2 (TelemacFile) second telemac-mascaret file
        @param record (int) Record to read data from
        @param var1 (str) Variable name of the first telemac-mascaret file
        @param var2 (str) Variable name of the second telemac-mascaret file
        @param var3 (str) Variable name of the third telemac-mascaret file. When
        not None, third and first var are multiplied
        @param var4 (str) Variable name of the fourth telemac-mascaret file. When
        not None, third and first var are multiplied
        @param relative (boolean) If true relative error is computed otherwise
                                  absolute one
        @param norm (str) Norm to apply (linf, l1, l2) default linf
        @param mass (np.array) mass matrix for errors computation
        @param check_name (bool) If False only check that we have the same
        number of variable but not necessaraly the same name
        """
        # If res1 or res2 is not specified, taking the same res file for both
        if res1 is None:
            res1 = res2
        if res2 is None:
            res2 = res1

        # Check if the two results have the same variables, records and points
        check_compatibility(res1, res2, record, check_name=check_name)

        # Computes errors:
        if var1 == 'all' or var2 == 'all':
            if var3 is not None:
                raise TelemacException(\
                        "Can't use var3 if var1 is not specified")
            if var4 is not None:
                raise TelemacException(\
                        "Can't use var4 if var2 is not specified")

            errs = []
            for idx, var in enumerate(res1.varnames):
                data1 = res1.get_data_value(var, record)
                data2 = res2.get_data_value(res2.varnames[idx], record)

                diff = compute_diff(data1, data2, relative=relative)
                err = compute_norm(diff, norm=norm, mass=mass)
                errs.append(err)
        else:
            if var1 not in res1.varnames:
                raise TelemacException("Variable 1 is not in result file")
            if var2 not in res2.varnames:
                raise TelemacException("Variable 2 is not in result file")

            data1 = res1.get_data_value(var1, record)
            data2 = res2.get_data_value(var2, record)

            if var3 is not None:
                if var3 not in res1.varnames:
                    raise TelemacException("Variable 1 is not in result file")
                data3 = res1.get_data_value(var3, record)
                data1 *= data3

            if var4 is not None:
                if var4 not in res1.varnames:
                    raise TelemacException("Variable 1 is not in result file")
                data4 = res2.get_data_value(var4, record)
                data2 *= data4

            diff = compute_diff(data1, data2, relative=relative)
            err = compute_norm(diff, norm=norm, mass=mass)
            errs = err

        # Return errors
        return errs

    def check_epsilons(self, str1, str2, eps, record=-1,
                       var1='all', var2='all',
                       relative=False, norm='linf', mass=None,
                       check_name=True):
        """
        Check the value of file1 - file2 for each variable or for given
        variables var1 and var2. Both result files must have the same mesh,
        the same variables and the same records. str1 and str2 can either be:
        - Name of a file
        - reference to a study name_study:T2DRES or name_study
        var1 and var2 can either be:
        - Name of a variable
        - 'all', all variables in file1 and file2 are checked

        @param str1 (str) Name of the first telemac-mascaret file
        @param str2 (str) Name of the second telemac-mascaret file
        @param eps (list) List of epsilon for each (if seq_and_par is True eps
        @param record (int) Record to read data from
        @param var1 (str) Variable name of the first telemac-mascaret file
        @param var2 (str) Variable of the second telemac-mascaret file
        @param relative (boolean) If true relative error is computed otherwise
                                  absolute one
        @param norm (str) Norm to apply (linf, l1, l2) default linf
        @param mass (np.array) mass matrix for errors computation
        @param check_name (bool) If False only check that we have the same
        number of variable but not necessaraly the same name
        """
        # TODO: Remove force sucess (was used ti bypass unset epsilons)
        # Converting string into file
        file1 = str1
        file2 = str2
        if ':' in str1:
            file1 = self.get_study_file(str1)
        if ':' in str2:
            file2 = self.get_study_file(str2)

        if file1 != file2:
            print(" "*8+"+> checking epsilon between files {} and {}:"\
                  .format(file1.replace(self.case_dir, '.'),
                          file2.replace(self.case_dir, '.')))
        else:
            print(" "*8+"+> checking epsilon for file {}:"\
                  .format(file1.replace(self.case_dir, '.')))

        # Opening the two telemac result files
        res1 = TelemacFile(file1)
        res2 = TelemacFile(file2)
        nvar = len(res1.varnames)

        # Compute errors:
        errs = self.compute_errors(res1=res1, res2=res2,
                                   var1=var1, var2=var2,
                                   record=record,
                                   relative=relative,
                                   norm=norm, mass=mass,
                                   check_name=check_name)
        # Epsilon check:
        force_success = False

        if var1 == 'all' or var2 == 'all':
            # Setting default espilon values
            # if eps is empty setting default value i.e. 1.e-4
            if eps == []:
                force_success = True
                print(" "*8+"~> TODO: Set Epsilon value")
                final_eps = [1.e-4]*nvar
            # If only one is given, it is used for all the variables
            elif len(eps) == 1:
                final_eps = [eps[0]] * nvar
            # Otherwise one must be given for each variable
            elif len(eps) != nvar:
                raise TelemacException(\
                    "Error in length of espilon is {} should be {}"\
                    .format(len(eps), nvar))
            else:
                final_eps = eps

            # Check epsilon for each variable in file1 and file2:
            for idx, eps in enumerate(final_eps):
                var = res1.varnames[idx]
                err = errs[idx]
                print(" "*10 + "- Difference for variable {}: {} (eps={})"\
                      .format(var, err, eps))
                passed = (err < eps)
                if not(passed or force_success):
                    raise TelemacException(\
                            "Epsilon reached in {} vs {}".format(file1, file2))
        else:
            # Setting default espilon values
            eps = eps[0]

            # Check epsilon
            err = errs
            print(" " * 10 + "- Difference for variables {}/{}: {} (eps={})"\
                  .format(var1, var2, err, eps))
            passed = (err < eps)
            if not(passed or force_success):
                raise TelemacException(\
                        "Epsilon reached in {} vs {}".format(file1, file2))
            errs = [err]

        # Closing the two files
        res1.close()
        res2.close()

        # Return errors
        return errs

    def get_ntimestep(self, str1):
        """
        Returns the number of time steps in a file
        """
        file1 = str1
        if ':' in str1:
            file1 = self.get_study_file(str1)

        res = TelemacFile(file1)
        ntimestep = res.ntimestep
        res.close()

        return ntimestep

    def pre(self):
        """
        Pre-treatment step
        """
        # Setting default values for action
        self.action_time['pre'] = [False, 0.0]
        start_time = time.time()
        # Running user defined pre-treatment
        self._pre()
        end_time = time.time()
        # Updating action_time information
        self.action_time['pre'] = [True, end_time - start_time]

    def run(self):
        """
        Run all the studies in the class
        """
        for name, study in self.studies.items():
            # Setting default values for action
            self.action_time[name] = [False, 0.0]
            start_time = time.time()
            self._run_case(name, study)
            end_time = time.time()
            # Updating action_time information
            self.action_time[name] = [True, end_time - start_time]

        for name, command in self.commands.items():
            # Running bash command
            self.action_time[name] = [False, 0.0]
            cfg = CFGS.configs[CFGS.cfgname]
            # If HPC config and not forced loca and command is hpc
            # Running command through slurm
            start_time = time.time()
            if cfg['HPC'] != {} and not self.options.mpi and command['hpc']:
                print("\n   ~> Submitting {}: \n{}\n".format(name, command['cmd']))
                hpc_cmd = get_hpc_cmd(cfg['HPC']).replace('<root>', cfg['root'])
                custom_path = 'cmd:{}:{}'.format(self.file_name, name)

                if 'id_log' in self.options:
                    hpc_cmd = hpc_cmd.replace('<id_log>', self.options.id_log
                                                          + ' ' + custom_path)
                else:
                    hpc_cmd = hpc_cmd.replace('<id_log>', 'id.log '+custom_path)
                stdin_file = cfg['HPC']['STDIN'][0]    # only one key for now
                stdin = cfg['HPC']['STDIN'][1]
                stdin = stdin.replace('<root>', cfg['root'])
                stdin = stdin.replace('<configName>', CFGS.cfgname)
                stdin = stdin.replace('<ncsize>', "1")
                stdin = stdin.replace('<nctile>', str(self.options.nctile))
                stdin = stdin.replace('<ncnode>', "1")
                stdin = stdin.replace('<email>', self.options.email)
                stdin = stdin.replace('<jobname>', name)
                ttime = time.strftime("%Y-%m-%d-%Hh%Mmin%Ss", time.localtime())
                stdin = stdin.replace('<time>', ttime)
                stdin = stdin.replace('<queue>', self.options.hpc_queue)
                stdin = stdin.replace('<walltime>', self.options.walltime)
                stdin = stdin.replace('\n ', '\n')
                # Setting command
                stdin = stdin.replace('<mpi_cmdexec>', command['cmd'])
                stdin = stdin.replace('<py_runcode>', command['cmd'])

                # ~~> Write to HPC_STDIN
                put_file_content(stdin_file, stdin.split('\n'))


                ret = subprocess.check_call(hpc_cmd, shell=True,
                                            executable=environ.get('SHELL', None))
                if ret != 0:
                    raise TelemacException(\
                        "The command below crashed: \n{}".format(hpc_cmd))
            else:
                print("\n   ~> Running {}: \n{}\n".format(name, command['cmd']))

                ret = subprocess.check_call(command['cmd'], shell=True,
                                            executable=environ.get('SHELL', None))

                if ret != 0:
                    raise TelemacException(\
                        "The command below crashed: \n{}".format(command['cmd']))
            end_time = time.time()
            # Updating action_time information
            self.action_time[name] = [True, end_time - start_time]

    def check_results(self):
        """
        Running verification and validation on the results of the studies
        """
        # Setting default values for action
        self.action_time['vnv'] = [False, 0.0]
        start_time = time.time()
        # Running user defined pre-treatment
        self._check_results()
        end_time = time.time()
        # Updating action_time information
        self.action_time['vnv'] = [True, end_time - start_time]

    def post(self):
        """
        Post treatment on the results of the studies
        """
        # Setting default values for action
        self.action_time['post'] = [False, 0.0]
        start_time = time.time()
        # Running user defined pre-treatment
        self._post()
        end_time = time.time()
        # Updating action_time information
        self.action_time['post'] = [True, end_time - start_time]

    @abstractmethod
    def _init(self):
        """
        Defines the general parameters (rank, tags)
        """
        pass

    @abstractmethod
    def _pre(self):
        """
        Defines the studies
        """
        pass

    @abstractmethod
    def _check_results(self):
        """
        Running comparaison and check between results file and reference file
        """
        pass

    @abstractmethod
    def _post(self):
        """
        Post treatment on the results of the studies
        """
        pass
