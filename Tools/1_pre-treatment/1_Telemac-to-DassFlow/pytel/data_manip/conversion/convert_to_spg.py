#!/usr/bin/env python
r"""@author Sebastien E. Bourban

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
        Creates a binary sponge file from a global model (SLF form)
          by interpolating on a given GEO model domain. The time series in
          the SPG file are extracted only at the nodes where the SPONGE mask
          value is more than 0.5.
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
__author__ = "Sebastien E. Bourban"
__date__ = "$02-Feb-2015 15:09:48$"


def main():
    """ Main function of convertToSPG """
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
Creates a binary sponge file from a global model (SLF form)
    by interpolating on a given GEO model domain. The time series in
    the SPG file are extracted only at the nodes where the SPONGE mask
    value is more than 0.5.
        '''),
        usage=' (--help for help)\n---------\n       =>  '
              '%(prog)s  geo-mask.slf in-result.slf out-sponge.slf\
               \n---------')
    parser.add_argument("args", default='', nargs=3)
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf new mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    geo_file = options.args[0]
    if not path.exists(geo_file):
        raise TelemacException(
                'Could not find geo_file: {}\n\n'.format(geo_file))

    # Read the new GEO file and its SPONGE mask variable
    print('   +> getting hold of the GEO file and of its SPONGE mask')
    geo = Selafin(geo_file)
    _, spg = geo.get_variables_at(0, subset_variables_slf(
              "BOTTOM: ;SPONGE mask: ", geo.varnames)[0])[0:2]
    print('   +> extracting the masked elements')
    # Keeping only masked nodes
    array_1d = np.in1d(geo.ikle2, np.sort(np.where(spg > 0.)[0]))
    mask = geo.ikle2[np.where(np.sum(array_1d.reshape(geo.nelem2, geo.ndp2),
                                     axis=1) > 0)]
    bor = np.unique(mask) + 1

    # Find corresponding (x,y) for the mask
    print('   +> getting hold of the GEO file and of its bathymetry')
    xys = np.vstack((geo.meshx[bor-1], geo.meshy[bor-1])).T

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf existing res ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    slf_file = options.args[1]
    if not path.exists(geo_file):
        raise TelemacException(
                'Could not find slf_file: {}\n\n'.format(slf_file))
    slf = Selafin(slf_file)
    slf.set_kd_tree()
    slf.set_mpl_tri()

    print('   +> support extraction')
    # Extract triangles and weigths in 2D
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
    support3d = zip(support2d, len(xys)*[range(slf.nplan)])

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes BND header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bnd_file = options.args[2]
    bnd = Selafin('')
    bnd.fole = {}
    bnd.fole.update({'hook': open(bnd_file, 'wb')})
    bnd.fole.update({'name': bnd_file})
    bnd.fole.update({'endian': ">"})     # big endian
    bnd.fole.update({'float': ('f', 4)})  # single precision

    # Meta data and variable names
    bnd.title = ''
    bnd.nbv1 = 5
    # /!\ ELEVATION has to be the first variable
    # (for possible vertical re-interpolation within TELEMAC)
    bnd.varnames = ['ELEVATION Z     ',
                    'VELOCITY U      ', 'VELOCITY V      ',
                    'SALINITY        ', 'TEMPERATURE     ']
    bnd.varunits = ['M               ',
                    'M/S             ', 'M/S             ',
                    '                ', '                ']
    bnd.nvar = bnd.nbv1
    bnd.varindex = range(bnd.nvar)

    # Sizes and mesh connectivity
    bnd.nplan = slf.nplan
    bnd.ndp2 = 3
    bnd.ndp3 = 6
    bnd.npoin2 = len(bor)
    bnd.npoin3 = bnd.npoin2*slf.nplan
    bnd.iparam = [0, 0, 0, 0, 0, 0, bnd.nplan, 0, 0, 0]
    bnd.ipob2 = bor   # /!\ Note that ipobO keeps the original numbering
    print('   +> masking and setting connectivity')
    # Set the array that only includes elements of geo.ikle2
    # with at least two nodes in bor
    ikle2 = mask
    # ~~> re-numbering ikle2 as a local connectivity matrix
    knolg, _ = np.unique(np.ravel(ikle2), return_index=True)
    knogl = dict(zip(knolg, range(len(knolg))))
    bnd.ikle2 = - np.ones_like(ikle2, dtype=np.int)
    for k in range(len(ikle2)):
        # /!\ bnd.ikle2 has a local numbering, fit to the boundary elements
        bnd.ikle2[k] = [knogl[ikle2[k][0]],
                        knogl[ikle2[k][1]],
                        knogl[ikle2[k][2]]]
    # Last few numbers
    bnd.nelem2 = len(bnd.ikle2)
    if slf.nplan > 1:
        bnd.nelem3 = bnd.nelem2*(slf.nplan-1)
    else:
        bnd.nelem3 = bnd.nelem2
    # 3D structures
    if slf.nplan > 1:
        bnd.ipob3 = np.ravel(np.add(np.repeat(bnd.ipob2, slf.nplan)
                                    .reshape((bnd.npoin2, slf.nplan)),
                                    bnd.npoin2*np.arange(slf.nplan)).T)
        bnd.ikle3 = \
            np.repeat(bnd.npoin2*np.arange(slf.nplan-1),
                      bnd.nelem2*bnd.ndp3)\
              .reshape((bnd.nelem2*(slf.nplan-1), bnd.ndp3)) + \
            np.tile(np.add(np.tile(bnd.ikle2, 2),
                           np.repeat(bnd.npoin2*np.arange(2), bnd.ndp2)),
                    (slf.nplan-1, 1))
    else:
        bnd.ipob3 = bnd.ipob2
        bnd.ikle3 = bnd.ikle2
    # Mesh coordinates
    bnd.meshx = geo.meshx[bor-1]
    bnd.meshy = geo.meshy[bor-1]

    print('   +> writing header')
    # Write header
    bnd.append_header_slf()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes BND core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('   +> setting variables')
    # TIME and DATE extraction
    bnd.tags['times'] = slf.tags['times']
    # VARIABLE extraction
    vrs = subset_variables_slf("ELEVATION Z: ;VELOCITY U: ;VELOCITY V: "
                               ";SALINITY: ;TEMPERATURE: ", slf.varnames)

    # Read / Write data, one time step at a time to support large files
    print('   +> reading / writing variables')
    pbar = ProgressBar(maxval=len(slf.tags['times'])).start()
    zeros = np.zeros((bnd.npoin3, 1), dtype=np.float)
    for time in range(len(slf.tags['times'])):
        data = get_value_history_slf(slf.file, slf.tags, [time], support3d,
                                     slf.nvar, slf.npoin3, slf.nplan, vrs)
        # special case for TEMPERATURE and SALINITY
        data[3] = np.maximum(data[3], zeros)
        data[4] = np.maximum(data[4], zeros)
        data = np.reshape(np.transpose(np.reshape(np.ravel(data),
                                                  (bnd.nvar, bnd.npoin2,
                                                   bnd.nplan)),
                                       (0, 2, 1)),
                          (bnd.nvar, bnd.npoin3))
        bnd.append_core_time_slf(time)
        bnd.append_core_vars_slf(data)
        pbar.update(time)
    pbar.finish()

    # Close bnd_file
    bnd.fole['hook'].close()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
