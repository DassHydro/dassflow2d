#!/usr/bin/env python3
"""@author TELEMAC-MASCARET CONSORTIUM

   @brief handles the compilation of the system
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
# ~~> dependencies towards other pytel/modules
from utils.messages import banner, svn_banner
from compilation.compil_tools import compile_api_files, \
    update_cmdf, compile_cmdf
from config import add_config_argument, update_config, CFGS

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#
DEBUG = False

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban; Noemie Durand"
__date__ = "$19-Jul-2010 08:51:29$"


def compute_config(rescan, bypass):
    """
    Update the configuration with the information needed by compileTELEMAC

    @param rescan (boolean) If true run rescan
    @param bypass (boolean) If True bypass error
    """
    CFGS.compute_compilation_info(\
                    rescan=rescan,
                    bypass=bypass)

def main():
    """ Main function of compileTELEMAC """

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n' + 72 * '~' + '\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
Compile the TELEMAC system:\n
1. rescan the tree dependencies if necessary
2. check which files need re-compilation
3. create object files, libraries, executable, and other binaries
    depending on your configuration settings
Work with all active configurations.
        '''))

    parser = add_config_argument(parser)
    parser.add_argument(
        "-m", "--modules", metavar="modules",
        dest="modules", default='',
        help="specify the list modules . separated, default is taken from config file")
    parser.add_argument(
        "-b", "--bypass", action="store_true",
        dest="bypass", default=False,
        help="will bypass execution failures and try to carry on "
             "(final report at the end)")
    parser.add_argument(
        "--rescan", action="store_true",
        dest="rescan", default=False,
        help="will redo the scan of sources for an update of "
             "all the cmdf files")
    parser.add_argument(
        "--clean", action="store_true",
        dest="cleanup", default=False,
        help="will erase all object, executable libraries from folder "
             "on the selected configs/modules")
    parser.add_argument(
        "-j", type=int,
        dest="ncsize", default=0,
        help="set the number of core used for the parallel "
             "compilation of objects")
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        dest="verbose", default=False,
        help="If given will print every command")
    options = parser.parse_args()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Environment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    update_config(options)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ banners ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    svn_banner(CFGS.get_root())

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Separately dealing with rescan for all configs ? ~~~~~~~~~~~~
    for cfgname in CFGS.configs:
        # Setting configuration name
        CFGS.cfgname = cfgname
        # Getting configuration
        cfg = CFGS.configs[cfgname]

        if options.cleanup:
            CFGS.clean_install(cfgname)

        print('\n\n' + '\n'.join(banner(cfgname)))
        print('Scanning the source code for:\n' + '~' * 72 + '\n')

        compute_config(options.cleanup and options.rescan,
                       options.bypass)
        CFGS.light_dump()
        # Only if we ask for a scan
        if options.rescan:
            update_cmdf(options.bypass, options.cleanup, options.verbose)

        # /!\ multiple configurations will now generate multiple rescan
        # (because of tags and adds, specific to some configurations)
        compute_config(False, options.bypass)

        # ~~ Scans all cmdf files found in all modules ~~~~~~~~~~~~~~~~~~~~~
        # Specifying what module to compile
        if options.modules == '':
            modules = []
        else:
            modules = options.modules.split(".")
        compile_cmdf(options.ncsize, modules, options.verbose)

        # Compiling api if asked for
        cfg_opt = cfg.get('options', [])
        if 'api' in cfg_opt:
            compile_api_files(silent=not options.verbose,
                              static='static' in cfg_opt,
                              hermes_only='hermes_only' in cfg_opt)

    print('\n\nMy work is done\n\n')
    sys.exit(0)


if __name__ == "__main__":
    main()
