#!/usr/bin/env python3
"""
Main script for 1D/2D coupling
"""

import json
import os
import sys
from argparse import ArgumentParser
from telapy.coupling.long_cpl_driver import LongCplDriver

def check_json(fichier):
    """
    Read json file
    @param fichier (str) path json file of mascaret
    @return dico (dict) paramater of json file
    """
    dico = {}
    try:
        with open(fichier) as json_data:
            dico = json.load(json_data)
    except IOError:
        raise IOError('\n {} FILE MISSING '
                      'ARRET DU PROGRAMME. \n\n'.format(fichier))
    return dico


def long_cpl_run(model_dir, keep_exec):
    """
    Entry point for the longitudinal coupling standalone launcher

    @param model_dir (str) path of the coupling case root
    @param keep_exec (bool) toggle exec dir cleaning. Default False = erase it
    """
    test_path = os.path.realpath(model_dir)
    list_dico = {}
    if os.path.isdir(test_path):
        if os.path.isfile(os.path.join(test_path, 'coupling_description.py')):
            sys.path.append(test_path)
            imp = __import__('coupling_description')
            coupling_def = imp.coupling_def
            config_run = imp.config_run
            for key in coupling_def['1D'].keys():
                list_dico['config_{}'.format(key)] = eval('imp.config_{}'.format(key))

            for key in coupling_def['2D'].keys():
                list_dico['config_{}'.format(key)] = eval('imp.config_{}'.format(key))

        else:

            file = os.path.join(test_path, 'ConfigRun.json')
            if os.path.exists(file):
                config_run = check_json(file)
            else:

                raise NameError('\nDICO NOT FOUND {}. '
                                'PROGRAM STOP. \n\n'.format('ConfigRun.json'))

            file = os.path.join(test_path, 'CouplingDef.json')
            if os.path.exists(file):
                coupling_def = check_json(file)
            else:

                raise NameError('\nDICO NOT FOUND {}. '
                                'PROGRAM STOP. \n\n'.format('CouplingDef.json'))

            # 2D model
            liste_tel = list(coupling_def["2D"].keys())
            tel_config = "config_{}".format(liste_tel[0])
            name_file = tel_config + ".json"
            path_tel_file = os.path.join(test_path, name_file)
            if os.path.exists(path_tel_file):
                list_dico[tel_config] = check_json(path_tel_file)
            else:
                raise NameError('\nDICO NOT FOUND {}. '
                                'PROGRAM STOP. \n\n'.format(name_file))

            # 1D model
            list_mas = list(coupling_def["1D"].keys())
            for i_b, name in enumerate(list_mas):
                name_dict = "config_{}".format(name)
                name_file = name_dict + '.json'
                path_mas_file = os.path.join(test_path, name_file)
                if os.path.exists(path_mas_file):
                    list_dico[name_dict] = check_json(path_mas_file)
                else:
                    raise IOError("Missing config  {}".format(path_mas_file))
                i_b += 1

        if not (os.path.isfile(os.path.join(test_path, 'CouplingDef.json')) or
                os.path.isfile(os.path.join(test_path, 'coupling_description.py'))):
            raise ValueError(test_path + ' Is not a test case directory.\nSTOP')
    else:
        raise ValueError(test_path + ' Is not a directory.\nSTOP')

    cl_pc = LongCplDriver(test_path, coupling_def, list_dico, keep_exec)
    cl_pc(config_run)

def launch_coupling(nb1d):
    """
    Wrapper of the coupled models
    @param nb1d (int) number of 1d model instances
    """
    from telapy.coupling.mascaret_cpl import mascaret_cpl
    from telapy.coupling.telemac2d_cpl import telemac2d_cpl
    from mpi4py import MPI

    cpl_comm = MPI.COMM_WORLD
    rank = cpl_comm.Get_rank()
    if rank < nb1d:
        mascaret_cpl()
    else:
        telemac2d_cpl()

def main():
    """
    Main function for 1d/2d cpl
    """
    parser = ArgumentParser()
    subparser = parser.add_subparsers(\
            help='converter command to do', dest='command')

    # Adding options for long_cpl
    long_parser = subparser.add_parser(\
            'long',
            help='Running longitudinal coupling')

    long_parser.add_argument(\
            '--src',
            dest="model_dir",
            help="Directory of the model where is found \n"
                 "coupling_def.py file by default it is '.'",
            type=str,
            default='.')

    long_parser.add_argument(\
            '--keep-exec',
            dest="keep_exec",
            action="store_true",
            help="Preserve the EXEC directory after the run",
            default=False)

    # Adding options for the code launcher
    launch_parser = subparser.add_parser('launcher', \
                                         help='Launching the coupled components')

    launch_parser.add_argument(\
                               '--n1d',
                               dest="ncsize1d",
                               help="Number of 1d models to start",
                               type=int,
                               default=1)

    args = parser.parse_args()

    if args.command == 'long':
        long_cpl_run(args.model_dir, args.keep_exec)
    elif args.command == 'launcher':
        launch_coupling(args.ncsize1d)
    else:
        args.print_help()

    sys.exit(0)

if __name__ == "__main__":
    main()
