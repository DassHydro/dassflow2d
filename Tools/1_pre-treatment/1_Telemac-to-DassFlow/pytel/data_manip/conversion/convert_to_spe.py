#!/usr/bin/env python
r"""@author Juliette C.E. Parisi and Sebastien E. Bourban

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
        Creates a binary liquid boundary file from a global model (SLF form)
          by interpolating on a given GEO model domain. The time series in
          the BND file are extracted only at liquid nodes as defined in the
          CONLIM file.

    @history 02/12/2013 -- Juliette C.E. Parisi
        Created draft script to write the binary liquid boundary file
            on the basis of the HYCOM global model results

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
from struct import unpack
import sys
from os import path
import numpy as np
# ~~> dependencies towards other pytel scripts
from utils.progressbar import ProgressBar
from utils.files import put_file_content
from utils.exceptions import TelemacException
from data_manip.conversion.convert_utm import to_lat_long

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
__author__ = "Juliette C.E. Parisi"
__date__ = "$02-Dec-2013 15:09:48$"


def main():
    """ Main function of convertToSPE """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Dependencies towards other modules ~~~~~~~~~~~~~~~~~~~~~~~~~~
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from data_manip.formats.selafin import Selafin
    from data_manip.formats.conlim import Conlim
    from pretel.meshes import xys_locate_mesh

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nInterpreting command line options\n'+'~'*72+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
A script to map spectral outter model results, stored as SELAFIN files,
 onto the
    spatially and time varying boundary of a spatially contained SELAFIN file
    of your choosing (your MESH).
        '''),
        usage=' (--help for help)\n---------\n       => '
              ' %(prog)s  open-bound.cli open-bound.slf in-outer-geo.slf '
              'in-outer-spec.slf out-bound.slf \n---------')
    parser.add_argument(
        "--ll2utm", dest="ll2utm", default=None,
        help="assume outer file is in lat-long and open-bound file in UTM")
    parser.add_argument("args", default='', nargs=5)
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ cli+slf new mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cli_file = options.args[0]
    if not path.exists(cli_file):
        raise TelemacException(
                '... the provided cli_file does not seem '
                'to exist: {}\n\n'.format(cli_file))
    geo_file = options.args[1]
    if not path.exists(geo_file):
        raise TelemacException(
                '... the provided geo_file does not seem to exist: '
                '{}\n\n'.format(geo_file))

    # Read the new CLI file to get boundary node numbers
    print('   +> getting hold of the CONLIM file and of its liquid boundaries')
    cli = Conlim(cli_file)
    # Keeping only open boundary nodes
    bor = np.extract(cli.bor['lih'] != 2, cli.bor['n'])

    # Find corresponding (x,y) in corresponding new mesh
    print('   +> getting hold of the GEO file and of its bathymetry')
    geo = Selafin(geo_file)
    if options.ll2utm is not None:
        zone = int(options.ll2utm)
        x, y = to_lat_long(geo.meshx[bor-1], geo.meshy[bor-1], zone)
    else:
        x = geo.meshx[bor-1]
        y = geo.meshy[bor-1]
    xys = np.vstack((x, y)).T

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf+spe existing res ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    slf_file = options.args[2]
    if not path.exists(slf_file):
        raise TelemacException(
                '... the provided slf_file does not seem to exist: '
                '{}\n\n'.format(slf_file))
    slf = Selafin(slf_file)
    slf.set_kd_tree()
    slf.set_mpl_tri()
    spe_file = options.args[3]
    if not path.exists(spe_file):
        raise TelemacException(
                '... the provided slf_file does not seem to exist: '
                '{}\n\n'.format(spe_file))
    spe = Selafin(spe_file)

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

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes BND header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bnd_file = options.args[4]
    bnd = Selafin('')
    bnd.fole = {}
    bnd.fole.update({'hook': open(bnd_file, 'wb')})
    bnd.fole.update({'name': bnd_file})
    bnd.fole.update({'endian': ">"})     # big endian
    bnd.fole.update({'float': ('f', 4)})  # single precision

    # Meta data and variable names
    bnd.title = spe.title
    # spectrum for new locations / nodes
    for i in range(len(bor)):
        bnd.varnames.append(('F'+('00'+str(i))[-2:]+' PT2D'+('000000' +
                             str(bor[i]))[-6:]+'                ')[:16])
        bnd.varunits.append('UI              ')
    bnd.nbv1 = len(bnd.varnames)
    bnd.nvar = bnd.nbv1
    bnd.varindex = range(bnd.nvar)

    # sizes and mesh connectivity / spectrum
    bnd.nplan = spe.nplan
    bnd.ndp2 = spe.ndp2
    bnd.ndp3 = bnd.ndp2
    bnd.npoin2 = spe.npoin2
    bnd.npoin3 = spe.npoin3
    bnd.iparam = spe.iparam
    bnd.ipob2 = spe.ipob2
    bnd.ikle2 = spe.ikle2
    # Last few numbers
    bnd.nelem2 = len(bnd.ikle2)
    bnd.nelem3 = bnd.nelem2
    bnd.ipob3 = bnd.ipob2
    bnd.ikle3 = bnd.ikle2
    # Mesh coordinates
    bnd.meshx = spe.meshx
    bnd.meshy = spe.meshy

    print('   +> writing header')
    # Write header
    bnd.append_header_slf()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes BND core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('   +> setting variables')
    # TIME and DATE extraction
    bnd.datetime = spe.datetime
    bnd.tags['times'] = spe.tags['times']

    # pointer initialisation
    f = spe.file['hook']
    endian = spe.file['endian']
    ftype, fsize = spe.file['float']

    # Identofy variables (required for support2d geo-locations)
    specloc = []
    for n, _ in support2d:
        specloc.extend(n)
    vars_indexes = np.unique(specloc)
    if fsize == 4:
        z = np.zeros((len(vars_indexes), spe.npoin2), dtype=np.float32)
        data = np.zeros(spe.npoin2, dtype=np.float32)
    else:
        z = np.zeros((len(vars_indexes), spe.npoin2), dtype=np.float64)
        data = np.zeros(spe.npoin2, dtype=np.float64)

    # Read / Write data, one time step at a time to support large files
    print('   +> reading / writing variables')
    pbar = ProgressBar(maxval=len(spe.tags['times'])).start()
    for itime in range(len(spe.tags['times'])):
        f.seek(spe.tags['cores'][itime])  # [itime] is the frame
        # to be extracted
        f.seek(4+fsize+4, 1)  # the file pointer is initialised
        bnd.append_core_time_slf(itime)

        # Extract relevant spectrum, where
        #  vars_indexes only contains the relevant nodes
        #  jvar varies from 0 to len(vars_indexes)
        jvar = 0
        for ivar in range(spe.nvar):
            # the file pointer advances through all records to keep on track
            f.seek(4, 1)
            if ivar in vars_indexes:
                z[jvar, :] = unpack(endian+str(spe.npoin2) +
                                    ftype, f.read(fsize*spe.npoin2))
                jvar += 1
            else:
                # the file pointer advances through
                # all records to keep on track
                f.seek(fsize*spe.npoin2, 1)
            f.seek(4, 1)

        # linear interpolation
        ivar = 0
        for b_n, l_n in support2d:
            data[:] = 0.
            for inod in range(len(b_n)):
                jvar = np.where(vars_indexes == b_n[inod])[0][0]
                data += l_n[inod]*z[jvar, :]
            bnd.append_core_vars_slf([data])
            ivar += 1

        pbar.update(itime)
    pbar.finish()

    # Close bnd_file
    bnd.fole['hook'].close()

    print('   +> writing out the file with coordinate to impose')
    dat = [str(len(bor)) + ' 0']
    for i in np.sort(bor):
        dat.append(str(i) + ' ' + repr(geo.meshx[i-1]) + ' ' +
                   repr(geo.meshy[i-1]) + ' 0.0')
    put_file_content(bnd_file+'.dat', dat)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
