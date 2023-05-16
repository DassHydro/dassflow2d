#!/usr/bin/env python3
r"""
    @author TELEMAC-MASCARET Consortium
    @brief Function to run a steering file using the api
"""
# Class Telemac2d import
from execution.telemac_cas import TelemacCas
from telapy.api.t2d import Telemac2d
from telapy.api.t3d import Telemac3d
from telapy.api.art import Artemis
from telapy.api.wac import Tomawac
from mpi4py import MPI
from argparse import ArgumentParser
import string
from os import path, environ


SCRIPT_TEMPLATE = """\
#!/usr/bin/env python3
# Class {module} import
from __future__ import print_function
from telapy.api.{short} import {module}
from mpi4py import MPI
# Creation of the instance Telemac2d
STEERING_FILE = "{steering_file}"
USER_FORTRAN = "{fortran_file}"
COMM = MPI.COMM_WORLD
# Initialising {module} instance
STUDY = {module}(STEERING_FILE, user_fortran=USER_FORTRAN, comm=COMM)
# Reading steering file
STUDY.set_case()
# Doing initialisation
STUDY.init_state_default()
# Time step loop
STUDY.run_all_time_steps()
# Ending computation
STUDY.finalize()
# Deleting {module} instance
del STUDY
"""

SHORT = {'telemac2d':'t2d',
         'telemac3d':'t3d',
         'tomawac':'wac',
         'sisyphe':'sis',
         'artemis':'art'}

def get_fortran_file(steering_file, module):
    """
    Get the fortran file from a cas (looks into coupled steering files as well)

    @param steering_file Name of the steering file
    @param module Name of the module
    """
    dico = path.join(environ['HOMETEL'], 'sources', module, module+'.dico')
    cas = TelemacCas(steering_file, dico)
    fortran_file = cas.get('FORTRAN FILE', '')

    if fortran_file == '':
        # Only searching in coupled files for telemac2d and telemac3d
        fortran_file = None
        if module in ["telemac2d", "telemac3d"]:
            cpl_with = cas.get('COUPLING WITH', '')
            if cpl_with == '':
                return None
            cpl_mods = cpl_with.lower().split(';')
            for cpl_mod in cpl_mods:
                cpl_dico = path.join(environ['HOMETEL'],
                                     'sources',
                                     cpl_mod,
                                     cpl_mod+'.dico')
                cpl_steering_file = cas.get(cpl_mod.upper()+' STEERING FILE')
                # Some coupled module do not have a dictionary (nestor, waqtel)
                if path.exists(cpl_dico):
                    cpl_cas = TelemacCas(cpl_steering_file, cpl_dico)
                    fortran_file = cpl_cas.get('FORTRAN FILE', '')
                    del cpl_cas
                    if fortran_file != '':
                        return fortran_file
            return None


    return fortran_file


def run(module, steering_file, stdout, log):
    """
    Running a full study

    @param module (str) Name of the module
    @param steering_file (str) Name of the steering file
    @param stdout (int) Output for TelApy (6 normal listing -1 into file)
    @param log (str) Logging level for TelApy (INFO, DEBUG)
    """
    # Creation of the instance Telemac2d
    comm = MPI.COMM_WORLD
    fortran = get_fortran_file(steering_file, module)
    if module == "telemac2d":
        study = Telemac2d(steering_file, user_fortran=fortran,
                          comm=comm, stdout=stdout, log_lvl=log)
    elif module == "telemac3d":
        study = Telemac3d(steering_file, user_fortran=fortran,
                          comm=comm, stdout=stdout, log_lvl=log)
    elif module == "artemis":
        study = Artemis(steering_file, user_fortran=fortran,
                        comm=comm, stdout=stdout, log_lvl=log)
    elif module == "tomawac":
        study = Tomawac(steering_file, user_fortran=fortran,
                        comm=comm, stdout=stdout, log_lvl=log)
    # Running telemac
    study.set_case()
    study.init_state_default()
    study.run_all_time_steps()
    comm.Barrier()
    study.finalize()
    # Instance delete
    del study

def dump_script(module, steering_file, fortran_file, script_file):
    """
    dump a api run into a file

    @param module (string) Name of the module
    @param steering_file (string) Name of the steering file
    @param fortran_file (string) Name of the fortran file
    @param script_file (string) Name of file in which we write the script
    """

    script = SCRIPT_TEMPLATE.format(\
            steering_file=steering_file,
            fortran_file=None if fortran_file == '' else fortran_file,
            module=string.capwords(module),
            short=SHORT[module])

    with open(script_file, 'w') as fobj:
        fobj.write(script)

if __name__ == "__main__":
    # Define a parser for the program options
    PARSER = ArgumentParser()
    PARSER.add_argument(\
             "module",
             choices=['telemac2d', 'telemac3d', 'artemis', 'tomawac'],
             help="name of the steering file")
    PARSER.add_argument(\
             "steering_file",
             help="name of the steering file")
    PARSER.add_argument(\
             "--double-run",
             dest="double_run",
             action="store_true",
             help="Running main computation twice")
    PARSER.add_argument(\
             "-v", "--verbose",
             dest="verbose",
             action="store_true",
             help="Display Telemac listing")
    PARSER.add_argument(\
             "--log",
             dest="log",
             default='INFO',
             choices=['INFO', 'DEBUG'],
             help="TelApy log level")
    PARSER.add_argument(\
             "-o", "--output-script",
             dest="output_script",
             default="",
             help="Will generate a python script running the case")
    # reading the options
    ARGS = PARSER.parse_args()

    if ARGS.output_script != '':
        dump_script(ARGS.module, ARGS.steering_file,
                    get_fortran_file(ARGS.steering_file, ARGS.module),
                    ARGS.output_script)
    else:
        STDOUT = 6 if ARGS.verbose else 0
        run(ARGS.module, ARGS.steering_file, STDOUT, ARGS.log)
        print("First run passed")
        if ARGS.double_run:
            run(ARGS.module, ARGS.steering_file, STDOUT, ARGS.log)
            print("Second run passed")

    print("My work is done")
