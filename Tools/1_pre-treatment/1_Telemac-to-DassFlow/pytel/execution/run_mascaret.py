#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium
"""
# ~~> dependencies towards standard python
import sys
import shutil
from os import path, remove
from argparse import ArgumentParser, RawDescriptionHelpFormatter
# ~~> dependencies towards other pytel/modules
from utils.messages import Messages, svn_banner
from config import add_config_argument, update_config, CFGS


def create_mascaret_files(cfg, cas):
    """
        Creates if not there the following files in the current folder:
        FichierCas.txt That contains the name of the steering file
        Abaques.txt ???
        Controle.txt ???

        param cfg Configuration object
        param cas Name of the cas file given as argument
    """
    # Always creating FichierCas.txt
    if path.isfile("FichierCas.txt"):
        remove("FichierCas.txt")

    print('~+> Creating FichierCas.txt')
    with open("FichierCas.txt", 'w') as fobj:
        fobj.write("'"+cas+"'\n")

    # Copying the abaque file if necessary
    if not path.isfile("Abaques.txt"):
        print('~+> Copying Abaques.txt')
        shutil.copyfile(path.join(cfg['root'], 'sources', 'mascaret', 'data',
                                  'Abaques.txt'),
                        "Abaques.txt")
    # Copying the controle file if necessary
    if not path.isfile("Controle.txt"):
        print('~+> Copying Controle.txt')
        shutil.copyfile(path.join(cfg['root'], 'sources', 'mascaret', 'data',
                                  'Controle.txt'),
                        "Controle.txt")
    # Copying Damocle Courlis dictionnary if necessary
    if not path.isfile("dico_Courlis.txt"):
        print('~+> Copying dico_Courlis.txt')
        shutil.copyfile(path.join(cfg['root'], 'sources', 'mascaret', 'data',
                                  'dico_Courlis.txt'),
                        "dico_Courlis.txt")


def run_mascaret():
    """
    Main function that runs the mascaret executable in the current folder
    """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+72*'~'+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n\
Run the mascaret executable in the current folder, given a CAS file.
        '''))
    parser.add_argument("args", nargs='*')
    # ~~> Environment
    parser = add_config_argument(parser)
    parser.add_argument(
        "-s", "--sortiefile", action="store_true", dest="sortie_file",
        default=False,
        help="specify whether there is a sortie file, default is no")
    parser.add_argument(
        "-b", "--bypass", action="store_true",
        dest="bypass", default=False,
        help="will bypass execution failures and try to carry on "
             "(final report at the end)")
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Environment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    update_config(options)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ banners ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    svn_banner(CFGS.get_root())

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for one configuration only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(options.args) < 1:
        print('\nThe name of the CAS file is required\n')
        parser.print_help()
        sys.exit(1)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads command line arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cas = options.args[0]
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for only one configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    xcpts = Messages()

    # still in lower case
    # parsing for proper naming
    CFGS.compute_execution_info()
    cfg = CFGS.configs[CFGS.cfgname]

    create_mascaret_files(cfg, cas)

    mascaret_exe = path.join(cfg['root'], 'builds', CFGS.cfgname, 'bin',
                             'mascaret' + cfg['sfx_exe'])
    _, code = xcpts.run_cmd(mascaret_exe, options.bypass)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reporting errors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if xcpts.not_empty() or code != 0:
        print('\n\nHummm ... I could not complete my work.\n'+'~'*72
              + xcpts.except_messages())
        sys.exit(1)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else:
        print('\n\nMy work is done\n\n')
        sys.exit(0)
