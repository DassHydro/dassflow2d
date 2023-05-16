#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium
   @brief Run the recollection step
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
import sys
from os import path
import argparse
# ~~> dependencies towards the root of pytel
# ~~> dependencies towards other pytel/modules
from execution.run import run_gretel
from execution.get import get_gretel_cmd
from config import add_config_argument, update_config, CFGS


def main():
    """ Main function of gretel """
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+72*'~'+'\n')
    parser = argparse.ArgumentParser(\
            description='Run the merging step (gretel)')
    parser = add_config_argument(parser)
    parser.add_argument(\
              "--input-file",
              dest="input_file",
              default='',
              help="Name of gretel parameter file (GRETEL.PAR). "\
                     "This option will surcharge all the others")
    parser.add_argument(\
              "--geo-file",
              dest="geo_file",
              default='T2DGEO',
              help="Name of the geometry file associated with the "\
                     "file to be merged")
    parser.add_argument(\
              "--geo-file-format",
              dest="geo_file_fmt",
              default='SERAFIN',
              help="Format of the geometry file(SERAFIN,SERAFIND or MED), "\
                     "default is SERAFIN")
    parser.add_argument(\
              "--res-file",
              dest="res_file",
              default='T2DRES',
              help="Name of the file to be merged")
    parser.add_argument(\
              "--res-file-format",
              dest="res_file_fmt",
              default='SERAFIN',
              help="Format of the geometry file(SERAFIN,SERAFIND or MED), "\
                     "default is SERAFIN")
    parser.add_argument(\
              "--bnd-file",
              dest="bnd_file",
              default='T2DCLI',
              help="Name of the boundary file")
    parser.add_argument(\
              "--ncsize",
              dest="ncsize",
              default=8,
              help="Number of partitions (should be equal to number of "\
                     "parallel processors), default is 8")
    parser.add_argument(\
              "--nplan",
              dest="nplan",
              default=0,
              help="Number of horizontal levels ,default is 0")
    parser.add_argument(\
              "--method",
              dest="method",
              default=1,
              help="Method for merging data ,default is 1")
    args = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Environment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    update_config(args)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for all configurations unless specified ~~~~~~~~~~~~~~~

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reporting errors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for all configurations unless specified ~~~~~~~~~~~~~~~
    CFGS.compute_execution_info()

    if args.input_file != "":
        with open(args.input_file, 'r') as f:
            geo_file = f.readline().strip('\n')
            geo_file_fmt = f.readline().strip('\n')
            bnd = f.readline().strip('\n')
            res_file = f.readline().strip('\n')
            res_file_fmt = f.readline().strip('\n')
            ncsize = f.readline().strip('\n')
            nplan = f.readline().strip('\n')
            method = f.readline().strip('\n')
    else:
        geo_file = args.geo_file
        geo_file_fmt = args.geo_file_fmt
        bnd = args.bnd_file
        res_file = args.res_file
        res_file_fmt = args.res_file_fmt
        ncsize = args.ncsize
        nplan = args.nplan
        method = args.method

    # Getting partel command from configuration
    pbin = path.join(CFGS.get_root(), 'builds', CFGS.cfgname, 'bin')
    grecmd = get_gretel_cmd(pbin, CFGS.configs[CFGS.cfgname])
    # Running paritionning

    run_gretel(grecmd, res_file, res_file_fmt, geo_file, geo_file_fmt, bnd,
               ncsize, nplan, method, False)

    print('\n\nMy work is done\n\n')
    sys.exit(0)

if __name__ == "__main__":
    main()
