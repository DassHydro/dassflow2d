#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

@brief Run the partiotionning step
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
import sys
from os import path
import argparse
# ~~> dependencies towards other pytel/modules
from execution.run import run_partel
from execution.get import get_partel_cmd
from config import add_config_argument, update_config, CFGS


def main():
    """ Main function of partel.py """
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+72*'~'+'\n')
    parser = argparse.ArgumentParser(
        description='Run the partitionning step (partel)')
    parser = add_config_argument(parser)
    parser.add_argument(
        "--input-file",
        dest="input_file",
        default='',
        help="Name of partel parameter file (PARTEL.PAR)."
             " This option will surcharge all the others")
    parser.add_argument(
        "--file",
        dest="geo_file",
        default='T2DGEO',
        help="Name of the file to be partitionned")
    parser.add_argument(
        "--file-format",
        dest="geo_file_fmt",
        default='SERAFIN',
        help="Format of the geometry file(SERAFIN,SERAFIND or MED), "
             "default is SERAFIN")
    parser.add_argument(
        "--bnd-file",
        dest="bnd_file",
        default='T2DCLI',
        help="Name of the boundary file associated to the mesh file, "
             "default is T2DCLI")
    parser.add_argument(
        "--ncsize",
        dest="ncsize",
        default=8,
        help="Number of partitions (should be equal to number of "
             "parallel processors), default is 8")
    parser.add_argument(
        "--section-name",
        dest="section_file",
        default='',
        help="Name of the section file, default no section file")
    parser.add_argument(
        "--zone-name",
        dest="zone_file",
        default='',
        help="Name of the zone file, default no zone file")
    parser.add_argument(
        "--weir-name",
        dest="weir_file",
        default='',
        help="Name of the weir file, default no weir file")
    parser.add_argument(
        "--partitioning-method",
        dest="part_method",
        default=1,
        help="Method used for the partitionning (1:metis, 2:scotch)")
    parser.add_argument(
        "--concat",
        dest="concat",
        action="store_true",
        default=False,
        help="If true concatenate partel output")
    parser.add_argument(
            "--mpi",
            dest="mpi",
            action="store_true",
            default=False,
            help="Run partel as executable (note using command given in systel.cfg)")
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
            bnd_file = f.readline().strip('\n')
            ncsize = f.readline().strip('\n')
            part_method = f.readline().strip('\n')
            section_file = f.readline().strip('\n')
            zone_file = f.readline().strip('\n')
            weir_file = f.readline().strip('\n')
            _ = f.readline().strip('\n')
            _ = f.readline().strip('\n')
            concat = f.readline().strip('\n')
    else:
        concat = 'YES' if args.concat else 'NO'
        geo_file = args.geo_file
        geo_file_fmt = args.geo_file_fmt
        bnd_file = args.bnd_file
        ncsize = args.ncsize
        section_file = args.section_file
        zone_file = args.zone_file
        weir_file = args.weir_file
        part_method = args.part_method

    # Getting partel command from configuration
    pbin = path.join(CFGS.get_root(), 'builds', CFGS.cfgname, 'bin')
    if args.mpi:
        exe_ext = CFGS.configs[CFGS.cfgname]['SYSTEM']['sfx_exe']
        parcmd = path.join(pbin, 'partel'+exe_ext+\
                                 ' < <partel.par> >> <partel.log>')
    else:
        parcmd = get_partel_cmd(pbin, CFGS.configs[CFGS.cfgname], '')
    # Running paritionning

    run_partel(parcmd, geo_file, geo_file_fmt, bnd_file, ncsize, False,
               section_file, zone_file, weir_file,
               geo_file, geo_file_fmt, part_method, concat)

    print('\n\nMy work is done\n\n')
    sys.exit(0)


if __name__ == "__main__":
    main()
