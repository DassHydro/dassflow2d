#!/usr/bin/env python3
"""@author TELEMAC-MASCARET Consortium
   @brief Refine a SERAFIN mesh
"""

import sys
import os
import subprocess as sp

CAS_REFINEMENT_CANVAS = \
"""
/
/ REFINEMENT OF MESH FILE USING STBTEL
/
/
/ INPUT FILE INFORMATION
/
UNIVERSAL FILE : '{inputFile}'
BOUNDARY UNIVERSAL FILE : '{inputBndFile}'
/
MESH GENERATOR : 'SELAFIN'
CUTTING ELEMENTS IN FOUR : YES
/
/ OUTPUT FILE INFORMATION
/
GEOMETRY FILE FOR TELEMAC : '{outputFile}'
BOUNDARY CONDITIONS FILE : '{outputBndFile}'
"""

def build_refine_cas(input_file, boundary_file, output_name):
    """
    Create the cas file for a refinement job

    @param input_file Name of the file to refine
    @param boundary_file Name of the boundary file associated with the file to
                         refine
    @param output_name Basename of the ouput

    @return (string) the steering case
    """
    return CAS_REFINEMENT_CANVAS.format(inputFile=input_file,
                                        inputBndFile=boundary_file,
                                        outputFile=output_name+".slf",
                                        outputBndFile=output_name+".cli")

def run_refine(input_file, output_file, root_dir, bnd_file):
    """
    Run a refinement using stbtel

    @param input_file (string) Name of the input file
    @param output_file (string) Name of the output_file
    @param root_dir (string) Path to the root of Telemac
    @param bnd_file (string) Boundary file
    """

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Identifying input and output informations ~~~~~~~~~~~~~~~~~~~


    # Treatment in case we are doing a refinement
    output_name, _ = os.path.splitext(output_file)
    cas = build_refine_cas(input_file, bnd_file, output_name)
    # Writting the steering file
    cas_name = 'stb.cas'
    with open(cas_name, "w") as fobj:
        fobj.write(cas)
    # Running stbtel
    path_stbtel = "stbtel.py"
    if root_dir is not None:
        path_stbtel = os.path.join(root_dir, "scripts",
                                   "python3", path_stbtel)
    stbtel_args = [path_stbtel, cas_name, "--mpi"]
    if root_dir is not None:
        stbtel_args += ["-r", root_dir]
    print("Calling: " + " ".join(stbtel_args))
    code = sp.call(stbtel_args)

    if code != 0:
        sys.exit(code)
    else:
        # Remove the case file
        os.remove(cas_name)

def stbtel_refine_parser(subparser):
    """
    Adding argument to parser for stbtel refinment

    @param subparser (ArgumentParser) the parser to update

    @return (ArgumentParser) the updated parser
    """

    parser = subparser.add_parser('refine',\
            help='Refinment of the mesh using stbtel')
    parser.add_argument(
        "input_file", default="",
        help="name of the input file also defines the input format")
    # output name option
    parser.add_argument(
        dest="output_file", default="",
        help="name of the output file also defines the output format")
    # the boundary file option
    parser.add_argument(
        "-b", "--boundary-file",
        dest="bnd_file", default="",
        help="name of the boundary file")
    # root directory
    parser.add_argument(
        "-r", "--root-dir",
        dest="root_dir", default=None,
        help="specify the root, default is taken from config file")

    return subparser
