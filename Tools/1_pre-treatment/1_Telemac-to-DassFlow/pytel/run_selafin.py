#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium

   Manipulation of selafin files
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from argparse import ArgumentParser
# ~~> dependencies towards other modules
from pretel.manip_selafin import scan, spec, chop, alter, merge, diff, \
                                 sample, subdivide, tesselate, calcs
from utils.exceptions import TelemacException

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
def add_arg(parser, name):
    """
    Add the argument name to parser

    @param parser (argparse.ArgumentParser) The parser
    @param name (str) name of the argument to add
    """
    if name == 'xvars':
        parser.add_argument(\
            "-v", "--vars",
            dest="xvars", default=None,
            help="specify which variables should remain (','-delimited)")
    elif name == 'core':
        parser.add_argument(\
            "-c", "--core", action="store_true",
            dest="core", default=False,
            help="specify whether to print statistics on the core "
                 "variables")
    elif name == 'time':
        parser.add_argument(\
            "-f", "--from",
            dest="tfrom", default="1",
            help="specify the first frame included")
        parser.add_argument(\
            "-s", "--stop",
            dest="tstop", default="-1",
            help="specify the last frame included "
                 "(negative from the end)")
        parser.add_argument(\
            "-d", "--step",
            dest="tstep", default="1",
            help="specify the step for the extraction of frames")
    elif name == 'replace':
        parser.add_argument(\
            "-r", "--replace", action="store_true",
            dest="freplace", default=False,
            help="if present, the output file will eventualy replace "\
                 "the input file")
    elif name == 'atitle':
        parser.add_argument(\
            "--title",
            dest="atitle", default=None,
            help="set the title of the SLF")
    elif name == 'areset':
        parser.add_argument(\
            "--reset", action="store_true",
            dest="areset", default=False,
            help="reset AT to zero second")
    elif name == 'adate':
        parser.add_argument(\
            "--date",
            dest="adate", default=None,
            help="set the start date of the SLF (dd-mm-yyyy)")
    elif name == 'atime':
        parser.add_argument(\
            "--time",
            dest="atime", default=None,
            help="set the start time of the SLF (hh:mm:ss)")
    elif name == 'eswitch':
        parser.add_argument(\
            "--endian", action="store_true",
            dest="eswitch", default=False,
            help="switch between endian encoddings")
    elif name == 'fswitch':
        parser.add_argument(\
            "--float", action="store_true",
            dest="fswitch", default=False,
            help="switch between DOUBLE and SINGLE precision float")
    elif name == 'aswitch':
        parser.add_argument(\
            "--switch", action="store_true",
            dest="aswitch", default=False,
            help="switch between VARIABLES and CLANDESTINES")
    elif name == 'aname':
        parser.add_argument(\
            "--name",
            dest="aname", default=None,
            help="change the name of a VARIABLE: 'OLD VAR=NEW VAR'")
    elif name == 'modif_time':
        parser.add_argument(\
            "--T+?",
            dest="atp", default="0",
            help="adds to the ATs")
        parser.add_argument(\
            "--T*?",
            dest="atm", default="1",
            help="scales the ATs")
    elif name == 'convert':
        parser.add_argument(\
            "--sph2ll",
            dest="sph2ll", default=None,
            help="convert from spherical to longitude-latitude")
        parser.add_argument(\
            "--ll2sph",
            dest="ll2sph", default=None,
            help="convert from longitude-latitude to spherical")
        parser.add_argument(\
            "--ll2utm",
            dest="ll2utm", default=None,
            help="zone info convert from longitude-latitude to UTM (giving XXX "
            "will let the script pick the zone otherwise specify zone with "
            "number + letter (for example 24S))")
        parser.add_argument(\
            "--utm2ll",
            dest="utm2ll", default=None,
            help="zone info to convert from UTM to longitude-latitude (for "
            "example 24S)")
    elif name == 'modif_coord':
        parser.add_argument(\
            "--X+?",
            dest="axp", default="0",
            help="adds to the meshx")
        parser.add_argument(\
            "--X*?",
            dest="axm", default="1",
            help="scales the meshx")
        parser.add_argument(\
            "--Y+?",
            dest="ayp", default="0",
            help="adds to the meshy")
        parser.add_argument(\
            "--Y*?",
            dest="aym", default="1",
            help="scales the meshy")
    elif name == 'modif_var':
        parser.add_argument(\
            "--Z?",
            dest="azname", default=None,
            help="will filter Z+ znd Z* ((Z* * val )+ Z+) "\
                 "operations on that VARIABLE name")
        parser.add_argument(\
            "--Z+?",
            dest="azp", default="0",
            help="adds to the VARIABLE")
        parser.add_argument(\
            "--Z*?",
            dest="azm", default="1",
            help="scales the VARIABLE")
    elif name == 'accuracy':
        parser.add_argument(\
            "--accuracy",
            dest="accuracy", default="5",
            help="significant figures for text display")
    elif name == 'points':
        parser.add_argument(\
            "--points",
            dest="points", default=None,
            help="extract data only at those locations")
        parser.add_argument(\
            "--nodes",
            dest="nodes", default=None,
            help="extract data only at those nodes")
    elif name == 'parallel':
        parser.add_argument(\
            "--parallel", action="store_true",
            dest="parallel", default=False,
            help="if option there, will assume input files have not been "\
                  "recollected, in which case you also need one example "\
                  "of the global file")
    else:
        print("Unknow options to add to parser: {}".format(name))
        print("Blame is on the developer")
        sys.exit(1)

    return parser

def scan_parser(subparser):
    """ Generate parser fo scan action """

    parser = subparser.add_parser('scan',\
        help='will print information about the SELAFIN, such as variables,'\
             'their vales etc.')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'core')
    parser = add_arg(parser, 'time')
    parser.add_argument("input_files", metavar='SELAFIN file', nargs="+")

    return subparser

def spec_parser(subparser):
    """ Generate parser fo spec action """
    parser = subparser.add_parser('spec',\
            help='will print information about a spectral file (also SELAFIN),'\
                 'such as frequencies, periodes, etc.')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'core')
    parser = add_arg(parser, 'time')
    parser = add_arg(parser, 'accuracy')
    parser.add_argument("input_files", metavar='Spectrum SELAFIN file',
                        nargs="+")

    return subparser

def chop_parser(subparser):
    """ Generate parser fo chop action """
    parser = subparser.add_parser('chop',\
       help='will chop a SELAFIN given a new set of time range and step (but'\
            ' alter is better)')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'time')
    parser = add_arg(parser, 'eswitch')
    parser = add_arg(parser, 'fswitch')
    parser = add_arg(parser, 'replace')
    parser = add_arg(parser, 'parallel')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser

def alter_parser(subparser):
    """ Generate parser fo alter action """
    parser = subparser.add_parser('alter',\
            help='will alter a SELAFIN file, choping or modifying time,'\
                 'converting its coordinates, extracting variables, etc')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'time')
    parser = add_arg(parser, 'replace')
    parser = add_arg(parser, 'parallel')
    parser = add_arg(parser, 'atitle')
    parser = add_arg(parser, 'areset')
    parser = add_arg(parser, 'adate')
    parser = add_arg(parser, 'atime')
    parser = add_arg(parser, 'aswitch')
    parser = add_arg(parser, 'eswitch')
    parser = add_arg(parser, 'fswitch')
    parser = add_arg(parser, 'aname')
    parser = add_arg(parser, 'modif_time')
    parser = add_arg(parser, 'modif_coord')
    parser = add_arg(parser, 'modif_var')
    parser = add_arg(parser, 'convert')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser

def merge_parser(subparser):
    """ Genrate parser for merge """
    parser = subparser.add_parser('merge',\
            help='will merge two files together, whether they are continuous'\
                 'simulations (same variables) or putting variables together'\
                 ' (same time definition)')
    parser = add_arg(parser, 'parallel')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'atitle')
    parser = add_arg(parser, 'areset')
    parser = add_arg(parser, 'adate')
    parser = add_arg(parser, 'atime')
    parser = add_arg(parser, 'aswitch')
    parser = add_arg(parser, 'eswitch')
    parser = add_arg(parser, 'fswitch')
    parser = add_arg(parser, 'aname')
    parser = add_arg(parser, 'modif_time')
    parser = add_arg(parser, 'modif_coord')
    parser = add_arg(parser, 'modif_var')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser


def diff_parser(subparser):
    """ Generate parser for diff """
    # TODO: set names for arguments
    parser = subparser.add_parser('diff',\
            help='Write the diff of all the variables of two files into '\
                 'a third one (the must have the same mesh)')
    parser.add_argument("args", metavar='SELAFIN file', nargs=3)

    return subparser

def sample_parser(subparser):
    """ Generate parser for sample """
    parser = subparser.add_parser('sample',\
            help='???')
    parser = add_arg(parser, 'parallel')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'time')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser

def subdivide_parser(subparser):
    """ Generate parser for subdivide """
    parser = subparser.add_parser('subdivide',\
            help='will subdivide a mesh by one iteration '\
            '(splitting all triangles in four others)'\
            ' deprecated use converter.py refine instead')
    parser = add_arg(parser, 'replace')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser

def tesselate_parser(subparser):
    """ Generate parser for tesselate """
    parser = subparser.add_parser('tesselate',\
            help='???')
    parser = add_arg(parser, 'replace')
    parser = add_arg(parser, 'convert')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser

def calcs_parser(subparser, name, help_msg):
    """ Generate parser for parser """
    parser = subparser.add_parser(name,\
            help=help_msg)
    parser = add_arg(parser, 'parallel')
    parser = add_arg(parser, 'xvars')
    parser = add_arg(parser, 'time')
    parser = add_arg(parser, 'modif_time')
    parser = add_arg(parser, 'modif_coord')
    parser = add_arg(parser, 'modif_var')
    parser = add_arg(parser, 'eswitch')
    parser = add_arg(parser, 'fswitch')
    parser.add_argument("args", metavar='SELAFIN file', nargs="+")

    return subparser

def main():
    """
    Options for each code name
        scan [--core] *.slf [*.slf]
        + '--core': print statistics on all variables, for each time step with
                      the core of all SELAFIN file present in args
        chop [--from] [--step] [--stop] in.slf out.slf
        + '--from': first frame included in out.slf
        + '--step': step used to extract the appropriate frame for out.slf
        + '--stop': last frame included in out.slf (unless step jumps over it)
          frame numbers (from and stop) being numbered from 0
        + '--vars': list those variables that are being extracted (all if empty)
        + '--replace': replace the input file by the output file, in which case
          multiple input files can be used
        alter [--title] [--date] [--time] ...
    """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+'~'*72+'\n')

    parser = ArgumentParser()
    subparser = parser.add_subparsers(\
            help='run_selafin commands to do', dest='command')

    subparser = chop_parser(subparser)
    subparser = scan_parser(subparser)
    subparser = spec_parser(subparser)
    subparser = alter_parser(subparser)
    subparser = merge_parser(subparser)
    subparser = diff_parser(subparser)
    subparser = calcs_parser(subparser, 'calcs', '???')
    subparser = calcs_parser(subparser, 'crunch', '???')
    subparser = calcs_parser(subparser, 'transf', '???')
    subparser = sample_parser(subparser)
    subparser = subdivide_parser(subparser)
    subparser = tesselate_parser(subparser)

    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads code name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if options.command == 'scan':
        scan(options)
    elif options.command == 'spec':
        spec(options)
    elif options.command == 'chop':
        chop(options)
    elif options.command == 'alter':
        alter(options)
    elif options.command == 'merge':
        merge(options)
    elif options.command == 'diff':
        diff(options)
    elif options.command == 'sample':
        sample(options)
    elif options.command in ['calcs', 'crunch', 'transf']:
        calcs(options, options.command)
    elif options.command == 'subdivide':
        subdivide(options)
    elif options.command == 'tessellate':
        tesselate(options)
    else:
        raise TelemacException(\
                '\nDo not know what to do with '
                'this code name: {}'.format(options.command))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)

if __name__ == "__main__":
    main()
