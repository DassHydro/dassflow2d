"""
Script for validation of the api
"""
from os import path, getcwd, chdir
from filecmp import cmp
import time
import shutil
from utils.exceptions import TelemacException
from execution.telemac_cas import TelemacCas, get_dico

def pre_api(my_vnv_study):
    """
    Duplicate study for api run

    @param my_vnv_study (vnv_study) The study in which to add api runs
    """

    _, old_time = my_vnv_study.action_time['pre']
    my_vnv_study.action_time['pre'] = [False, 0.0]
    start_time = time.time()

    for name, study in my_vnv_study.studies.items():
        api_name = name + "_api"

        # We need to copy all the files from study into the new study
        # Build directory for api study
        api_vnv_working_dir = my_vnv_study.build_vnv_working_dir(api_name)

        # No api for those modules
        if study.code_name in ['stbtel', 'postel3d']:
            print('  ~> Not doing api run for this module {}'
                  .format(study.code_name))

        # Temporary study just to copy the files
        study.copy_files(api_vnv_working_dir,
                         verbose=my_vnv_study.options.verbose,
                         copy_cas_file=True)

        cmd = "cd {wdir} && mpirun -n {ncsize} template.py {module} {cas} --double-run"\
                .format(wdir=api_vnv_working_dir,
                        ncsize=max(1, study.ncsize),
                        module=study.code_name,
                        cas=path.basename(study.steering_file))

        # Handle in coupling cases when the coupled steering file is using the
        # same file as the main one
        # This does not work with the api as we do not have a temporary folder
        if study.code_name != 'mascaret':
            cas = study.cas
            in_files = [cas.get(key) for key in cas.in_files]

            pwd = getcwd()

            chdir(api_vnv_working_dir)

            modified = False
            for mod, tmp_cpl_cas in study.cpl_cases.items():
                # Looking for input file that is in both the main steering and the
                # coupled one
                cpl_cas = TelemacCas(path.basename(tmp_cpl_cas.file_name),
                                     get_dico(mod))
                for key in cpl_cas.in_files:
                    ffile = cpl_cas.get(key)
                    if ffile in in_files:
                        root, ext = path.splitext(ffile)
                        new_file = root + '_' + mod +ext
                        print(" ~> Copying {} -> {}".format(ffile, new_file))
                        # Copying file
                        shutil.copy2(path.join(api_vnv_working_dir, ffile),
                                     path.join(api_vnv_working_dir, new_file))
                        print(" ~> Modifying in {}: {}".format(mod, key))
                        # Changing value in steering file
                        cpl_cas.set(key, new_file)

                        modified = True
                # If we made some modification overwritting the steering case
                if modified:
                    cpl_cas.write(path.join(api_vnv_working_dir, cpl_cas.file_name))
                del cpl_cas

            chdir(pwd)

        my_vnv_study.add_command(api_name, cmd, hpc=True)

    end_time = time.time()
    # Updating action_time information
    my_vnv_study.action_time['pre'] = [True, old_time + end_time - start_time]

def check_api(my_vnv_study):
    """
    Compare output files between normal run and api run

    @param my_vnv_study (vnv_study) The study for which to do the check
    """

    _, old_time = my_vnv_study.action_time['vnv']
    my_vnv_study.action_time['vnv'] = [False, 0.0]
    start_time = time.time()

    for name, study in my_vnv_study.studies.items():

        # No api for those modules
        if study.code_name in ['stbtel', 'postel3d']:
            print('  ~> Not doing api run for this module {}'
                  .format(study.code_name))

        api_name = name + "_api"

        api_vnv_working_dir = my_vnv_study.build_vnv_working_dir(api_name)
        vnv_working_dir = my_vnv_study.build_vnv_working_dir(name)

        # Building list of cas for which to check output files
        cases = [study.cas]

        for cas in study.cpl_cases.values():
            cases.append(cas)

        for cas in cases:
            for key in cas.out_files:
                out_file = cas.get(key)
                print(" "*8+"+> checking files for {} vs {}"\
                        .format(name, api_name))
                ffile = path.join(vnv_working_dir, out_file)
                api_ffile = path.join(api_vnv_working_dir, out_file)
                compare = cmp(ffile, api_ffile)

                if not compare:
                    raise TelemacException("Difference in {} ({} vs {})"\
                            .format(out_file, name, api_name))


    end_time = time.time()
    # Updating action_time information
    my_vnv_study.action_time['vnv'] = [True, old_time + end_time - start_time]
