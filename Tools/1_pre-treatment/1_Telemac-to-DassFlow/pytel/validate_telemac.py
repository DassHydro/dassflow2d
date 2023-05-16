#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium
   @brief
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
import time
from os import path, walk, chdir, remove, sep, listdir, getcwd
from argparse import ArgumentParser
from shutil import rmtree
# ~~> dependencies towards other pytel/modules
from utils.messages import Messages, svn_banner
from utils.files import check_sym_link
from utils.exceptions import TelemacException
from vvytel.run_notebook import run_notebook
from vvytel.vnv_api import check_api, pre_api
from vvytel.report_class import Report
from config import update_config, CFGS
from runcode import add_runcode_argument
from execution.wait_for_jobs import check_job_slurm, get_job_time_slurm

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

TAGS = ["telemac2d",
        "telemac3d",
        "mascaret",
        "courlis",
        "sisyphe",
        "nestor",
        "tomawac",
        "aed",
        "waqtel",
        "python2",
        "python3",
        "apistudy",
        "api_mascaret",
        "coupling",
        "postel3d",
        "stbtel",
        "khione",
        "artemis",
        "gaia",
        "gotm",
        "full_valid",
        "med",
        "fv",
        "api"]


def check_python_rank_tags(py_file, options):
    """
    Checks if a Python vnv script match the rank and tags options
    """
    val_dir = path.dirname(py_file)
    abs_py_file = path.abspath(py_file)

    chdir(val_dir)

    # Importing vnv_class from py_file
    try:
        # Code foor Python 3.5+
        import importlib.util
        # This allow Python script decalared in the example folder to be loaded
        sys.path.append(val_dir)
        spec = importlib.util.spec_from_file_location("vnv_module", py_file)
        vnv_stuff = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(vnv_stuff)
    except:
        from importlib.machinery import SourceFileLoader

        vnv_stuff = SourceFileLoader("vnv_module", py_file).load_module()


    name = path.splitext(py_file)[0]
    my_vnv_study = vnv_stuff.VnvStudy(name, abs_py_file, options)

    # Checkcing ranks if will run all ranks <= options.rank
    rank = my_vnv_study.rank
    rank_ok = rank <= options.rank

    if not rank_ok and options.verbose:
        print('\n ~> '+py_file)
        print('     > nothing to do here (rank):')
        print('       {} > {}: '.format(rank, options.rank))

    # Checking tags will run all the test contains one of the tags in
    # options.tags And skip the ones with a - before them (-mascaret will
    # skip file with mascaret tag)
    tags = my_vnv_study.tags
    if tags != []:
        tags_ok = False
        opt_tags = options.tags.split(',')

        if '+' in options.tags:
            for opt_tag in opt_tags:
                if '+' in opt_tag:
                    # All the tag in opt_tag must be in tag
                    tag_ok = True
                    for tag2 in opt_tag.split('+'):
                        tag_ok = tag_ok and (tag2 in tags)
                elif '-' in opt_tag:
                    # - means reverse tag must no be in
                    tag_ok = opt_tag[1:] not in tags
                else:
                    tag_ok = opt_tag in tags
                tags_ok = tags_ok or tag_ok
        else:
            for tag in tags:
                # If -tag that means that the Python should not be run
                if '-'+tag in opt_tags:
                    tags_ok = False
                    break
                tag_ok = tag in opt_tags
                # Checking that at least one of the tags is in opt_tags
                tags_ok = tags_ok or tag_ok
    else:
        raise TelemacException("Missing tag in Python file:\n"+py_file)

    if not tags_ok and options.verbose:
        print('\n ~> '+py_file)
        print('     > nothing to do here (tag):')
        print('       File tags: {}'.format(','.join(tags)))
        print('       Options tags: {}'.format(options.tags))

    # Cleaning up sys.path
    sys.path.remove(val_dir)

    return tags_ok and rank_ok

def get_list_python_files(cfg, options):
    """
    Get the list of files to run.
    Checking list given, rank, tags...

    @param cfg (Dict) Configuration information
    @param options (ArgumentParser) List of arguments

    @returns list_files (list) List of the files to run (full path)
    """

    list_files = []
    if options.args != []:
        for ffile in options.args:
            list_files.append(path.realpath(ffile))
    else:
        # Looping on all folders to find the scripts to execute
        # Loop on modules
        for code_name in cfg['VALIDATION']:
            val_root = cfg['val_root']

            dirpath, dirnames, _ = next(walk(path.join(val_root, code_name)))
            for ddir in dirnames:
                _, _, filenames = next(walk(path.join(dirpath, ddir)))
                for fle in filenames:
                    root, ext = path.splitext(path.basename(fle))
                    if ext == '.py' and root[0:4] == 'vnv_':
                        # check rank and tag
                        file_name = path.join(dirpath, ddir, fle)
                        if check_python_rank_tags(file_name, options):
                            list_files.append(path.realpath(file_name))

    return list_files

def run_validation_python_mpi(cfg, options, report, xcpts):
    """
    Run validation for vnv Python scripts Normale mode

    @param cfg (Dict) Configuration information
    @param options (ArgumentParser) List of arguments
    @param report (Report) Time of actions
    @param xcpts () Error handler
    """
    # Building list of files to run
    list_files = get_list_python_files(cfg, options)

    n_files = len(list_files)
    root = CFGS.get_root()
    for ifile, py_file in enumerate(sorted(list_files)):
        if options.cleanup or options.full_cleanup:
            clean_vnv_working_dir(py_file, full=options.full_cleanup)
        else:
            print('\n\nValidation < {}/{} > of {}'\
                  .format(ifile+1, n_files, py_file.replace(root, '<root>')))
            run_python(py_file, options, report, xcpts)

def run_validation_python_slurm(cfg, options, report, xcpts):
    """
    Run validation for vnv Python scripts slurm mode

    @param cfg (Dict) Configuration information
    @param options (ArgumentParser) List of arguments
    @param report (Report) Time of actions
    @param xcpts () Error handler
    """
    # Building list of files to run
    list_files = get_list_python_files(cfg, options)

    if list_files == []:
        print("Nothing to run (check tags and rank)")
        return

    n_files = len(list_files)
    root = CFGS.get_root()

    if options.cleanup or options.full_cleanup:
        for ifile, py_file in enumerate(sorted(list_files)):
            clean_vnv_working_dir(py_file, full=options.full_cleanup)

        return

    # Making sure that the file is not there before first run
    jobid_file = options.id_log
    if path.exists(jobid_file):
        remove(jobid_file)

    # ~> First submission run
    options.vnv_pre = True
    options.vnv_run = True
    options.vnv_post = False
    options.vnv_check = False
    options.bypass = True

    if options.hpc_queue == '':
        raise TelemacException(
            "Option --queue is mandatary with --vnv-mode=slurm")

    print("  ~> Submission part")
    for ifile, py_file in enumerate(sorted(list_files)):
        print('\n\nValidation < {}/{} > of {}'\
              .format(ifile+1, n_files, py_file.replace(root, '<root>')))
        run_python(py_file, options, report, xcpts)

    # Removing from list files all the ones that crashed in the first
    # submission run
    run_list_files = list_files.copy()
    for error in xcpts.messages:
        if error['name'] in run_list_files:
            run_list_files.remove(error['name'])

    # ~> Waiting for jobs to finish

    jobs = {}
    jobs_ini = {}
    run_times = {}
    crashed = {'failed':[], 'timeout':[]}

    # File is generated by the first run
    # In case no run was launched in the previous command
    if not path.exists(jobid_file):
        actual_len = 0
        start_time = time.time()
    else:
        # Building dictionary of jobs:
        with open(jobid_file, 'r') as f:
            for line in f.readlines():
                job_id, action_path = line.split(';')
                if job_id == '':
                    raise TelemacException(\
                       "Error in the job id file. "\
                       "Generated by hpc_runcode in systel.cfg:\n{}"
                       .format(jobid_file))
                tmp2 = action_path.replace('\n', '')
                if 'cmd:' in tmp2:
                    _, abs_py_file, action = tmp2.split(':')
                    working_dir = path.dirname(abs_py_file)
                else:
                    tmp = tmp2.split(sep)
                    if 'mascaret' in tmp or 'courlis' in tmp:
                        abs_py_file = sep.join(tmp[:-1])+".py"
                        action = tmp[-1]
                        working_dir = tmp2
                    else:
                        abs_py_file = sep.join(tmp[:-2])+".py"
                        action = tmp[-2]
                        working_dir = tmp2

                jobs[job_id] = {'name':abs_py_file+':'+action,
                                'py_file':abs_py_file,
                                'action':action,
                                'working_dir':working_dir}

        jobs_ini.update(jobs)

        # Check job status
        print("  ~> Waiting for completion")
        prev_len = 0
        actual_len = len(jobs)
        # Waiting time between each check in second
        wait_time = 10

        start_time = time.time()
        time.sleep(60)

    while actual_len != 0:
        # Only printing remaining jobs if there was some changes
        if prev_len != actual_len:
            print("Remaining jobs: ", len(jobs))
        t_start = time.time()
        for job_id in list(jobs.keys()):
            state = check_job_slurm(job_id)
            # Job crashed
            if state == 'failed':
                crashed['failed'].append(jobs[job_id])
                del jobs[job_id]
            # job timed out
            elif state == 'timeout':
                crashed['timeout'].append(jobs[job_id])
                del jobs[job_id]
            # Job is done
            elif state == 'success':
                run_time = get_job_time_slurm(job_id)
                run_times[jobs[job_id]['name']] = run_time
                del jobs[job_id]
            # Otherwise job is still running
        t_end = time.time()
        # Only wait if the loop was done in less than wait_time
        if (t_end - t_start) < wait_time:
            time.sleep(wait_time)
        # Update info on len
        prev_len = actual_len
        actual_len = len(jobs)

    elapsed_time = time.time() - start_time
    time_str = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    print("Waited {} for jobs to complete".format(time_str))

    # Adding run times to the report
    for action_name, run_time in run_times.items():
        abs_py_file, action = action_name.split(':')

        rank = report.values[abs_py_file]['pre']['rank']

        report.add_action(abs_py_file, rank,
                          action, run_time, True)

    # Building new list of files (without the ones that crashed)
    new_list_files = []
    for py_file in run_list_files:
        # Extract folder of validation from script name (minus estension)
        py_folder, _ = path.splitext(py_file)
        failed_jobs = crashed['timeout'] + crashed['failed']

        # Chekc if the file is in the failed files
        failed = False
        for job in failed_jobs:
            if job['py_file'] == py_file:
                failed = True
                break

        # If it is next file
        if failed:
            continue

        new_list_files.append(py_file)

    # Adding exception for all the run that crashed
    for crash_type, failed in crashed.items():
        if failed != []:
            for fail in failed:
                xcpts.add_messages([{'name':fail['name'],
                                     'msg':'The job {}'.format(crash_type)}])

    print("  ~> Displaying listing of all runs")
    # Displaying listings (before merging because merging will remove temporary folder)
    for val in jobs_ini.values():
        run_dir = val['working_dir']
        print('\n\nListing for {}:'\
              .format(val['name'].replace(path.realpath(root), '<root>')))
        # If cmdexec hpc mode listing is in the temporary folder
        if 'hpc_cmdexec' in cfg:
            for ffile in listdir(run_dir):
                if ffile[:4] == 'tmp_' and \
                   path.isdir(path.join(run_dir, ffile)):
                    run_dir = path.join(run_dir, ffile)
                    break

        for ffile in listdir(run_dir):
            if ffile.endswith(".out"):
                with open(path.join(run_dir, ffile), 'r',
                          encoding='utf-8') as f:
                    print(f.read())

    # If we are in hpc_cmdexec configuration (only out_telemac is in the batch
    # job)
    # Running on more pass to do the merge step
    # Second run
    if 'hpc_cmdexec' in cfg:
        print("  ~> Merging part")
        options.vnv_pre = True
        options.vnv_run = True
        options.vnv_post = False
        options.vnv_check = False
        options.bypass = True
        options.merge = True
        options.split = False
        options.run = False
        options.compileonly = False
        # Running only on jobs that finished
        for ifile, py_file in enumerate(sorted(new_list_files)):
            print('\n\nValidation < {}/{} > of {}'\
                  .format(ifile+1, n_files, py_file.replace(root, '<root>')))
            run_python(py_file, options, report, xcpts)



    # Second run
    options.vnv_pre = True
    options.vnv_run = False
    options.vnv_post = True
    options.vnv_check = True
    options.bypass = True
    print("  ~> Check + Post-traitment")
    # Running only on jobs that finished
    for ifile, py_file in enumerate(sorted(new_list_files)):
        print('\n\nValidation < {}/{} > of {}'\
              .format(ifile+1, n_files, py_file.replace(root, '<root>')))
        run_python(py_file, options, report, xcpts, time_from_report=True)


def clean_vnv_working_dir(py_file, full=False):
    """
    Clean the working directory

    @param py_file (str) Path of the vnv_study file
    @param full (boolean) If True just delete the entire folder
    """

    vnv_study_dir, _ = path.splitext(py_file)
    if full:
        rmtree(vnv_study_dir)
    else:
        for dirpath, dirnames, _ in walk(vnv_study_dir):
            for dirname in dirnames:
                if dirname == CFGS.cfgname:
                    rmtree(path.join(dirpath, dirname))


def run_python(py_file, options, report, xcpts, time_from_report=False):
    """
    Run a vnv Python script

    @param py_file (str) Name of the Python file to run
    @param options (ArgumentParser) Options of the script
    @param report (Report) Contains execution time
    @param xcpts () Error Handler
    @param time_from_report (bool) If true update vnv_study time array by what
    is in the report before running post
    """
    try:
        abs_py_file = path.abspath(py_file)
        if not path.isfile(abs_py_file):
            raise TelemacException(\
               '\nNot able to find your Python file:\n{}'\
               .format(abs_py_file))

        val_dir = path.dirname(abs_py_file)

        chdir(val_dir)
        # Importing vnv_class from py_file
        try:
            # Code foor Python 3.5+
            import importlib.util
            # This allow Python script decalared in the example folder to be
            # loaded
            sys.path.append(val_dir)
            spec = importlib.util.spec_from_file_location("vnv_module", py_file)
            vnv_stuff = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(vnv_stuff)
        except:
            from importlib.machinery import SourceFileLoader

            vnv_stuff = SourceFileLoader("vnv_module", py_file).load_module()

        name = path.splitext(py_file)[0]
        my_vnv_study = vnv_stuff.VnvStudy(name, abs_py_file, options)

        # Pre-treatment part
        # It is always done
        my_vnv_study.pre()

        if options.api:
            pre_api(my_vnv_study)

        # Cleaning ?
        if options.cleanup or options.full_cleanup:
            my_vnv_study.clean_vnv_working_dir(full=options.full_cleanup)
            return

        # Execution part
        if options.vnv_run:
            chdir(val_dir)
            my_vnv_study.run()

            # cleaning temporary files (if no post):
            if not options.vnv_post:
                for ffile in my_vnv_study.temporary_files:
                    remove(path.join(val_dir, ffile))

        # Check part
        if options.vnv_check:
            chdir(val_dir)
            my_vnv_study.check_results()
            if options.api:
                check_api(my_vnv_study)

        # Post_treatment part
        if options.vnv_post:
            # Update time from what if in the report
            if time_from_report:
                for action, val in report.values[abs_py_file].items():
                    my_vnv_study.action_time[action] = \
                            [val['passed'], val['time']]
            chdir(val_dir)
            my_vnv_study.post()

            # cleaning temporary files:
            for ffile in my_vnv_study.temporary_files:
                remove(path.join(val_dir, ffile))

    except Exception as exc:
        if options.bypass:
            xcpts.add_messages([{'name':py_file,
                                 'msg':str(exc)}])
        else:
            raise exc
    finally:
        #Updating report information
        if my_vnv_study is not None:
            for action, actions in my_vnv_study.action_time.items():
                report.add_action(abs_py_file, my_vnv_study.rank,
                                  action, actions[1], actions[0])
    # Cleaning up sys.path
    if val_dir in sys.path:
        sys.path.remove(val_dir)


def run_validation_notebooks(options, report, xcpts):
    """
    Run validation of the notebooks

    @param options (ArgumentParser) Options of the script
    @param report (Report) Contains execution time
    @param xcpts (Message) Error handler
    """
    root = CFGS.get_root()

    if options.args != []:
        # Files to run given in arguments
        nb_files = options.args
    else:
        # Looking through notebook folder for notebook to run
        nb_files = []

        for dirpath, _, ffiles in walk(path.join(root, 'notebooks')):
            for ffile in ffiles:
                # Skipping jupyter temporary folders
                if '.ipynb' in dirpath:
                    continue
                # If we have a notebook
                if '.ipynb' in ffile:
                    nb_files.append(path.join(dirpath, ffile))

    # Removing exlcuded notebooks
    if options.nb_exclude is not None:
        options.nb_exclude = options.nb_exclude.strip("'")
        for nb_file in list(nb_files):
            for exclude in options.nb_exclude.split(','):
                if nb_file.endswith(exclude+".ipynb"):
                    if nb_file in nb_files:
                        print("  ~> Excluding: ", nb_file)
                        nb_files.remove(nb_file)

    # Run notebook validation
    n_nb = len(nb_files)
    for i, nb_file in enumerate(sorted(nb_files)):
        print('Validation <{}/{}> ~> Running notebook {} in {}'\
              .format(i+1, n_nb, path.basename(nb_file), path.dirname(nb_file)))
        try:
            start = time.time()
            run_notebook(nb_file, options.nb_timeout,
                         update_nb=options.nb_update)
            end = time.time()
            report.add_notebook(nb_file, end-start, True)
        except Exception as exc:
            if options.bypass:
                report.add_notebook(nb_file, 0.0, False)
                xcpts.add_messages([{'name':nb_file,
                                     'msg':str(exc)}])
            else:
                raise exc

def set_parser():
    """
    Defining parser for validateTELEMAC.py

    @return (Namespace) Arguments values
    """

    print('\n\nLoading Options and Configurations\n' + 72 * '~' + '\n')
    parser = ArgumentParser( \
        description=('''\n
Used to validate the TELEMAC system against a benchmark of test cases for
a certain rank, and a certain tag'''))

    parser = add_runcode_argument(parser)
    parser.add_argument( \
        "-b", "--bypass", action="store_true", dest="bypass", default=False,
        help="will bypass execution failures and try to carry on "\
              "(final report at the end)")
    # Combine with all filters above, "rank" now controls everything
    # and Jenkins can control "rank"
    parser.add_argument( \
        "-k", "--rank", dest="rank", type=int, default=4,
        help="specify the ranks to be validated all rank lower or equal to "
             "the value will be run")
    parser.add_argument( \
        "--tags", dest="tags", default='all',
        help=\
         "specify tags (, separated) to run "\
         "  '-tag' will do the opposite and "\
         "tag1+tag2 will run cases that has both tag1 and tag2), "\
         "default is all of them")
    parser.add_argument( \
        "--valrootdir", dest="val_root", default='',
        help="specify the directory in which to search the validation cases, "\
              "default is taken from config file")
    ### vnv options
    parser.add_argument( \
        "--vnv-mode", dest="vnv_mode", default='mpi',
        choices=['mpi', 'slurm'],
        help=\
""" Defines validation execution mode (does not work for notebooks):
- mpi: Normal sequential run of all the validation Python files.
- slurm: Will run validation in 'cluster' mode.
         This means that first valdiation will be run with
           --vnv-run to submit all the jobs to the scheduler.
         Then we will wait for them to be finished.
         Then we run the second validation with --vnv-check and --vnv-post.""")
    parser.add_argument(
        "--id-log", dest="id_log", default=None,
        help="File containing id of submitted jobs (--vnv-mode=slurm)")
    parser.add_argument( \
        "--vnv-pre", action="store_true", dest="vnv_pre", default=False,
        help="Only do pre-treatment")
    parser.add_argument( \
        "--vnv-run", action="store_true", dest="vnv_run", default=False,
        help="Only do execution for each study")
    parser.add_argument( \
        "--vnv-check", action="store_true", dest="vnv_check", default=False,
        help="Only do check of results (epsilons)")
    parser.add_argument( \
        "--vnv-post", action="store_true", dest="vnv_post", default=False,
        help="Only do post-treatment")
    parser.add_argument( \
        "--report-name", dest="report_name", default='',
        help="will create a csv containing information on the validation "\
             "such as execution time, rank, if it passed...")
    parser.add_argument( \
        "--clean", action="store_true", dest="cleanup", default=False,
        help="will erase all object, executable, result files "\
              "from subfolders for the actual configuration")
    parser.add_argument( \
        "--full-clean", action="store_true", dest="full_cleanup", default=False,
        help="will erase all vnv study folders regarding of configurations")

    # Options for notebook
    parser.add_argument(
        "--notebook",
        dest="notebook",
        action="store_true", default=False,
        help="Run validation of notebook")
    parser.add_argument(
        "--notebook-timeout",
        dest="nb_timeout", type=int, default=60000,
        help="Time after whihc the notebook will be killed if still running")
    parser.add_argument(
        "--notebook-update",
        dest="nb_update",
        action="store_true", default=False,
        help="Update notebook file with the runned one")
    parser.add_argument(
        "--notebook-exclude",
        dest="nb_exclude",
        default=None,
        help=",separated list of notebook (name without extension for example index,) to exclude")


    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true", default=False,
        help="More verbose validation")

    # Options for api
    parser.add_argument(
        "--api",
        dest="api",
        action="store_true", default=False,
        help="Run validation of api")

    parser.add_argument("args", metavar='Python file(s)', nargs='*')
    options = parser.parse_args()

    if options.id_log is None:
        options.id_log = path.join(getcwd(), 'id.log')

    # Conversion of options.tags (replacing all by list) and checking that the
    # value is valid
    # Removing quotes
    tmp_tag = options.tags.strip("'\"").replace(';', ',')
    options.tags = tmp_tag
    # Checking that tags are valid
    for tag in options.tags.split(','):
        if '+' in tag:
            for and_tag in tag.split('+'):
                # Removing - if in tag
                ttag = and_tag[1:] if and_tag[0] == '-' else and_tag
                if ttag not in TAGS:
                    raise TelemacException(\
                       "Unknow tag: {tag}\nTags available: {tags}"\
                       .format(tag=ttag, tags=','.join(TAGS)))
        else:
            if tag == 'all':
                continue
            # Removing - if in tag
            ttag = tag[1:] if tag[0] == '-' else tag
            if ttag not in TAGS:
                raise TelemacException(\
                   "Unknow tag: {tag}\nTags available: {tags}"\
                   .format(tag=ttag, tags=','.join(TAGS)))

    # Replacing all by list of tags
    if 'all' in options.tags.split(','):
        options.tags = options.tags.replace('all', ','.join(TAGS))

    # If pre, run, post are all false switching them to true
    if not(options.vnv_pre or options.vnv_run or
           options.vnv_check or options.vnv_post):
        options.vnv_pre = True
        options.vnv_run = True
        options.vnv_check = True
        options.vnv_post = True
    else:
        # Options not available in this mode
        if options.vnv_mode == 'slurm':
            raise TelemacException(
                "option vnv-pre, vnv-run, vnv-check, "
                "vnv-post not available in this vnv-mode")

    return options

def config_corrections(options, cfgname):
    """
    overwrite configuration wiht options arguments

    @param options (Value) List of argparse options
    @param cfgname (string) Name of the configuration

    @return (dict) The configuration info
    """
    CFGS.cfgname = cfgname
    cfg = CFGS.configs[cfgname]
    # still in lower case
    if options.val_root != '':
        cfg['val_root'] = options.val_root
    # parsing for proper naming
    CFGS.compute_vnv_info()
    print('    +> ' + cfgname + ': ' + ', '.join(cfg['VALIDATION'].keys()))
    if options.cleanup:
        cfg['REBUILD'] = 2

    return cfg

def main():
    """
     @brief Main function of validateTELEMAC.
    """
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~ Handles input arguments
    options = set_parser()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Environment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    update_config(options)
    version = CFGS.configs[CFGS.cfgname].get('version', 'trunk')

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ banners ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    svn_banner(CFGS.get_root(), version)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Works for all configurations unless specified ~~~~~~~~~~~~~~~
    # Checking if symlink is available
    if options.use_link and not check_sym_link(options.use_link):
        raise TelemacException(\
                '\nThe symlink option is only available on Linux systems. '
                'Remove the option and try again')

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Forces not to use any Xwindows backend for Jenkins ~~~~~~~~~~
    if options.vnv_post:
        import matplotlib.pyplot as plt

        plt.switch_backend('agg')

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Reporting errors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xcpts = Messages()

    # ~~~~ Reporting summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if options.notebook:
        type_valid = 'notebooks'
    else:
        type_valid = 'examples'

    report = Report(options.report_name, type_valid)

    # ~~~ Running validation
    cfg = config_corrections(options, CFGS.cfgname)

    if options.notebook:
        run_validation_notebooks(options, report, xcpts)
    else:
        if options.vnv_mode == 'slurm':
            run_validation_python_slurm(cfg, options, report, xcpts)
        else:
            run_validation_python_mpi(cfg, options, report, xcpts)

    # Writting report
    if options.report_name != '':
        report.write()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Reporting errors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if xcpts.not_empty():
        print('\n\nHummm ... I could not complete my work.\n'
              '{}\n{}'.format('~' * 72, xcpts.except_messages()))
        sys.exit(1)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else:
        print('\n\nMy work is done\n\n')
        sys.exit(0)

if __name__ == "__main__":
    main()
