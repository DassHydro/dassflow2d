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
from config import add_config_argument, update_config, CFGS
from execution.telemac_cas import TelemacCas, get_dico
from utils.exceptions import TelemacException
from pretel.check_cas import check_cas
from pretel.diff_cas import diff_cas

def translate(cas):
    """
    Translating cas_file in french and english

    @param cas (TelemacCas) steering file
    """
    print("\n ~> Translating {}".format(path.basename(cas.file_name)))
    cas.write_fr_gb()


def sort(cas, keep_comments):
    """
    Sorting steering file using rubriques order

    @param cas (TelemacCas) steering file

    @parma keep_comments (bool) If true sort will append comments from original
    file at the end of the file
    """

    root, ext = path.splitext(cas.file_name)

    out_file = root + "_sorted" + ext

    print("\n ~> Sorting {}".format(path.basename(cas.file_name)))
    cas.write(out_file, keep_comments=keep_comments)



def add_module_arg(parser):
    """ adding module options """
    parser.add_argument(
        "module",
        choices=['postel3d', 'telemac2d', 'telemac3d', 'tomawac', 'artemis',
                 'sisyphe', 'gaia', 'waqtel', 'khione', 'stbtel'],
        help="Name of the module")
    return parser

def translate_parser(subparser):
    """ Generate parser for translate """
    parser = subparser.add_parser('translate',\
        help="Generate a french and english version of the steering file "
             "(cas_file suffixed with _fr and _gb)")
    parser = add_module_arg(parser)
    parser = add_config_argument(parser)
    parser.add_argument(
        "cas_file",
        help="Name of the steering file to translate")

    return subparser

def sort_parser(subparser):
    """ Generate parser for sort """
    parser = subparser.add_parser('sort',\
        help="Rewrites the steering file using rubriques to sort the keywords "
             "cas_file suffixed with _sorted")
    parser = add_module_arg(parser)
    parser = add_config_argument(parser)
    parser.add_argument(
        "cas_file",
        help="Name of the steering file to sort")

    parser.add_argument(
        "--keep-comments", action="store_true",
        dest="keep_comments", default=False,
        help="When sorting will append all original comments "
             "at the end of the file")

    return subparser

def diff_parser(subparser):
    """ Generate parser for diff """
    parser = subparser.add_parser('diff',\
        help="Diff between all steering files give as argument into a csv file."
             " By default only the different keyword will be in the file add "
             "option --all to have them all")
    parser = add_module_arg(parser)
    parser = add_config_argument(parser)
    parser.add_argument(
        "cas_files", nargs='+',
        help="Name of the steering file to compare")
    parser.add_argument(
        dest="csv_name",
        help="Name of the output csv file")

    parser.add_argument(
        "--delimiter",
        dest="delimiter", default=";",
        help="Delimiter in the csv file")

    parser.add_argument(
        "--all", action="store_true",
        dest="allkey", default=False,
        help="If given will When sorting will append all original comments "
             "at the end of the file")

    return subparser

def check_parser(subparser):
    """ Generate parser for check """
    parser = subparser.add_parser('check',\
        help="Check that content of steering file is valid (files exist, time "
             "are consistent between computation and input files)")
    parser = add_module_arg(parser)
    parser = add_config_argument(parser)
    parser.add_argument(
        "cas_file",
        help="Name of the steering file to check")

    return subparser


def main():
    """ Main function of manip_cas.py """
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+72*'~'+'\n')
    parser = argparse.ArgumentParser(
        description='Manipulation on steering files')
    subparser = parser.add_subparsers(\
        help='manip_cas commands to do', dest='command')
    subparser = translate_parser(subparser)
    subparser = sort_parser(subparser)
    subparser = diff_parser(subparser)
    subparser = check_parser(subparser)

    args = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Environment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    update_config(args)
    cfg = CFGS.configs[CFGS.cfgname]
    CFGS.compute_execution_info()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Searching for the dictionary associated with the steering case
    dico_file = get_dico(args.module)
    if not path.exists(dico_file):
        raise TelemacException(
            'Could not find the dictionary file: {}'.format(dico_file))

    if args.command == 'check':
        cas = TelemacCas(args.cas_file, dico_file, check_files=False)
        check_cas(args.module, cas)
    elif args.command == 'translate':
        cas = TelemacCas(args.cas_file, dico_file, check_files=False)
        translate(cas)
    elif args.command == 'sort':
        cas = TelemacCas(args.cas_file, dico_file, check_files=False)
        sort(cas, args.keep_comments)
    elif args.command == 'diff':
        diff_cas(args.cas_files, args.csv_name, dico_file,
                 delimiter=args.delimiter, allkey=args.allkey)
    else:
        print(cas)

    print('\n\nMy work is done\n\n')
    sys.exit(0)


if __name__ == "__main__":
    main()
