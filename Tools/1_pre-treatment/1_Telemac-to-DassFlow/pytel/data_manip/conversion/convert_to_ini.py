#!/usr/bin/env python
r"""@author Juliette C.E. Parisi, Michael S. Turnbull and Sebastien E. Bourban

   @note ... this work is based on a collaborative effort between
  .________.                                                          ,--.
  |        |                                                      .  (  (
  |,-.    /   HR Wallingford                EDF - LNHE           / \_ \_/ .--.
  /   \  /    Howbery Park,                 6, quai Watier       \   )   /_   )
   ,.  `'     Wallingford, Oxfordshire      78401 Cedex           `-'_  __ `--
  /  \   /    OX10 8BA, United Kingdom      Chatou, France        __/ \ \ `.
 /    `-'|    www.hrwallingford.com         innovation.edf.com   |    )  )  )
!________!                                                        `--'   `--

   @brief
      Beware the the variables number / names are set in the code below.
      Creates an initial condition file from a global model (SLF form)
        by interpolating on a given GEO model domain.

   @history 02/12/2013 -- Juliette C.E. Parisi
      Created draft script to write the binary liquid boundary file
         on the basis of the HYCOM global model results

   @history 14/02/2014 -- Michael S. Turnbull
      CFurther adaption and correction for application to HYCOM

   @history 11/11/2014 -- Sebastien E. Bourban
      Heavy modifictaions to make it generic and in order to
         correctly fill-in the ipobO and the ikle in both 2D and 3D
         based on fancy numpy programing.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
from os import path
import numpy as np
# ~~> dependencies towards other pytel scripts
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
__author__ = "Juliette C.E. Parisi"
__date__ = "$02-Dec-2013 15:09:48$"


def main():
    """ Main function of convertToIni """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Dependencies towards other modules ~~~~~~~~~~~~~~~~~~~~~~~~~~
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from data_manip.formats.selafin import Selafin
    from data_manip.extraction.parser_selafin import subset_variables_slf, \
        get_value_history_slf
    from pretel.meshes import xys_locate_mesh

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+'~'*72+'\n')
    parser = ArgumentParser(
       formatter_class=RawDescriptionHelpFormatter,
       description=('''\n
A script to map 2D or 3D outter model results stored in a SELAFIN file,\
 onto the
   one frame of contained SELAFIN file of your choosing (your MESH).
      '''),
       usage=' (--help for help)\n---------\n       =>  '
             '%(prog)s  geo-mesh.slf in-result.slf out-init.slf \n---------')
    parser.add_argument("args", default='', nargs=3)
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf new mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    geo_file = options.args[0]
    if not path.exists(geo_file):
        raise TelemacException(
            '... the provided geo_file does not seem to exist: {}'
            '\n\n'.format(geo_file))

# Find corresponding (x,y) in corresponding new mesh
    print('   +> getting hold of the GEO file and of its bathymetry')
    geo = Selafin(geo_file)
    xys = np.vstack((geo.meshx, geo.meshy)).T
    _ = geo.get_variables_at(0,
                             subset_variables_slf("BOTTOM: ",
                                                  geo.varnames)[0])[0]

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf existing res ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    slf_file = options.args[1]
    if not path.exists(slf_file):
        raise TelemacException(
            '... the provided geo_file does not seem to exist: {}'
            '\n\n'.format(slf_file))

    slf = Selafin(slf_file)
    slf.set_kd_tree()
    slf.set_mpl_tri()

    print('   +> support extraction')
    # Extract triangles and weights in 2D
    support2d = []
    ibar = 0
    pbar = ProgressBar(maxval=len(xys)).start()
    for xyi in xys:
        support2d.append(xys_locate_mesh(xyi, slf.ikle2, slf.meshx, slf.meshy,
                                         slf.tree, slf.neighbours))
        ibar += 1
        pbar.update(ibar)
    pbar.finish()
    # Extract support in 3D
    support3d = list(zip(support2d, len(xys)*[range(slf.nplan)]))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes INI header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ini_file = options.args[2]
    ini = Selafin('')
    ini.fole = {}
    ini.fole.update({'hook': open(ini_file, 'wb')})
    ini.fole.update({'name': ini_file})
    ini.fole.update({'endian': ">"})     # big endian
    ini.fole.update({'float': ('f', 4)})  # single precision

    # Meta data and variable names
    ini.title = ''
    ini.nbv1 = 5
    # /!\ ELEVATION has to be the first variable
    # (for possible vertical re-interpolation within TELEMAC)
    ini.varnames = ['ELEVATION Z     ',
                    'VELOCITY U      ', 'VELOCITY V      ',
                    'SALINITY        ', 'TEMPERATURE     ']
    ini.varunits = ['M               ',
                    'M/S             ', 'M/S             ',
                    '                ', '                ']
    ini.nvar = ini.nbv1
    ini.varindex = range(ini.nvar)

    # sizes and mesh connectivity
    ini.nplan = slf.nplan
    ini.ndp2 = 3
    ini.ndp3 = 6
    ini.npoin2 = geo.npoin2
    ini.npoin3 = geo.npoin2*ini.nplan
    ini.nelem2 = geo.nelem2
    ini.nelem3 = ini.nelem2*(ini.nplan-1)

    print('   +> setting connectivity')
    ini.ikle3 = \
        np.repeat(geo.npoin2*np.arange(ini.nplan-1),
                  geo.nelem2*ini.ndp3)\
        .reshape((geo.nelem2*(ini.nplan-1), ini.ndp3)) + \
        np.tile(np.add(np.tile(geo.ikle2, 2),
                np.repeat(geo.npoin2*np.arange(2), geo.ndp2)),
                (ini.nplan-1, 1))
    ini.ipob3 = np.ravel(np.add(np.repeat(geo.ipob2, ini.nplan)
                                  .reshape((geo.npoin2, ini.nplan)),
                                geo.npoin2*np.arange(ini.nplan)).T)
    ini.iparam = [0, 0, 0, 0, 0, 0, ini.nplan, 0, 0, 0]

    # Mesh coordinates
    ini.meshx = geo.meshx
    ini.meshy = geo.meshy

    print('   +> writing header')
    # Write header
    ini.append_header_slf()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes INI core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('   +> setting variables')
    ini.tags['times'] = slf.tags['times']
    # VARIABLE extraction
    vrs = subset_variables_slf("ELEVATION Z: ;VELOCITY U: ;VELOCITY V: "
                               ";SALINITY: ;TEMPERATURE: ", slf.varnames)

    # Read / Write data for first time step
    zeros = np.zeros((ini.npoin3, 1), dtype=np.float)

    print('   +> extracting variables')
    data = get_value_history_slf(slf.file, slf.tags, [0], support3d, slf.nvar,
                                 slf.npoin3, slf.nplan, vrs)

    # special case for TEMPERATURE and SALINITY
    data[3] = np.maximum(data[3], zeros)
    data[4] = np.maximum(data[4], zeros)
    print('   +> correcting variables')
    # duplicate values below bottom
    data = np.reshape(
            np.transpose(
              np.reshape(
                np.ravel(data), (ini.nvar, ini.npoin2, ini.nplan)),
              (0, 2, 1)),
            (ini.nvar, ini.npoin3))
    print('   +> writing variables')
    ini.append_core_time_slf(0)
    ini.append_core_vars_slf(data)

    # Close ini_file
    ini.fole['hook'].close()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
