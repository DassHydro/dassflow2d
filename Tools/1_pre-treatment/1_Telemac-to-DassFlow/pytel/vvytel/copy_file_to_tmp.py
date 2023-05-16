r"""@author TELEMAC-MASCARET Consortium

    @note ... this work is based on a collaborative effort between
                 TELEMAC-MASCARET consortium members
    @brief Function for validation of the Python API
"""

from os import chdir, mkdir, path, sep
from glob import glob
import shutil
from execution.telemac_cas import TelemacCas
###############################################################################
def copy_file_to_tmp(test_dir, tmp_dir, module, root_dir, skip_test):
    #@TODO: Merge with vvytel/copy_file_to_valid_dir
    """
       Copy all the files needed by the test case into the temporary folder.

       @param test_dir	path to the test case to validate
       @param tmp_dir	path to the test case temporary folder
       @param module	Name of the module
       @param root_dir	Root directory of the installation
       @param skip_test	Test cases to skip
    """
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if not path.exists(tmp_dir):
        mkdir(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        mkdir(tmp_dir)
    chdir(tmp_dir)
    # Getting list on input/output files from the dictionary
    dico_file = path.join(root_dir, 'sources', module, module+'.dico')
    # Getting list of steering file
    cas_files_path = glob(test_dir +sep+'*.cas')
    list_file = []
    for cas_file in cas_files_path:
        if path.basename(cas_file) in skip_test:
            continue
        shutil.copyfile(cas_file, path.basename(cas_file))
        cas = TelemacCas(cas_file, dico_file)
        user_fortran = None
        # Looping on input files
        for key in cas.in_files:
            ffile = cas.values[key]
            # if we have a user fortran
            if 'FORTRAN' in key:
                if path.exists(path.join(tmp_dir, ffile)) and\
                    not path.isfile(path.join(tmp_dir, ffile)):
                    shutil.rmtree(path.join(tmp_dir, ffile))
                user_fortran = ffile
                if path.isfile(test_dir+sep+ffile):
                    shutil.copyfile(test_dir+sep+ffile, ffile)
                else:
                    shutil.copytree(test_dir+sep+ffile, ffile)
            else:
                shutil.copyfile(test_dir+sep+ffile, ffile)
        list_file.append((path.basename(cas_file), user_fortran))
    return list_file

