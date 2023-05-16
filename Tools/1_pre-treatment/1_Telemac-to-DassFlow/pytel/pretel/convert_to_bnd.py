#!/usr/bin/env python
r"""@author TELEMAC-MASCARET Consortium

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
from os import path
import numpy as np
# ~~> dependencies towards other pytel scripts
from utils.progressbar import ProgressBar
from utils.exceptions import TelemacException
from data_manip.formats.selafin import Selafin
from data_manip.formats.conlim import Conlim
from data_manip.extraction.parser_selafin import subset_variables_slf, get_value_history_slf
from pretel.meshes import xys_locate_mesh

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
__author__ = "Juliette C.E. Parisi"
__date__ = "$02-Dec-2013 15:09:48$"

def generate_bnd_parser(subparser):
    """
    Add argument for generate_bnd to parser

    @param subparser (ArgumentParser) The paser to update

    @return (ArgumentParser) The updated parser
    """
    parser = subparser.add_parser('generate_bnd',\
            help='A script to map 2D or 3D outter model results, stored as'\
                 'SELAFIN files, onto the spatially and time varying boundary '\
                 'of a spatially contained SELAFIN file of your choosing (your'\
                 ' MESH)')
    parser.add_argument(\
        "cli_file", default='',
        help='Open boundary')
    parser.add_argument(\
        "geo_file", default='',
        help='Geometry file')
    parser.add_argument(\
        "slf_file", default='',
        help='Result file')
    parser.add_argument(\
        "bnd_file", default='',
        help='Ouput BND file')
    parser.add_argument("--varnames",\
        default='ELEVATION Z;VELOCITY U;VELOCITY V;SALINITY;TEMPERATURE',\
        help='list of variables to read from slf_file ; separated')
    parser.add_argument("--varunits",\
        default='M;M/S;M/S;;',\
        help='list of variables to read from slf_file ; separated')

    return subparser


def generate_bnd(cli_file, geo_file, slf_file, bnd_file, varnames, varunits,
                 showbar=True):
    """
    @param cli_file
    @param geo_file
    @param slf_file
    @param bnd_file
    @param varnames
    @param varunits
    @param showbar (boolean) If True display a showbar for the progress
    """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ cli+slf new mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not path.exists(cli_file):
        raise TelemacException(\
             '... the provided cli_file does not seem to exist:'
             ' {}\n\n'.format(cli_file))
    if not path.exists(geo_file):
        raise TelemacException(\
            '... the provided geo_file does not seem to exist: '
            '{}\n\n'.format(geo_file))

    if len(varnames) != len(varunits):
        raise TelemacException(\
          'Not the same number of variables and units\nvarnames: {}\nvarunits: {}'
          '{}\n\n'.format(varnames, varunits))


    # Read the new CLI file to get boundary node numbers
    print('   +> getting hold of the Conlim file and of its liquid boundaries')
    cli = Conlim(cli_file)
    # Keeping only open boundary nodes
    bor = np.extract(cli.bor['lih'] != 2, cli.bor['n'])

    # Find corresponding (x,y) in corresponding new mesh
    print('   +> getting hold of the GEO file and of its bathymetry')
    geo = Selafin(geo_file)
    xys = np.vstack((geo.meshx[bor-1], geo.meshy[bor-1])).T
    _ = geo.get_variables_at(0,\
                  subset_variables_slf("BOTTOM: ", geo.varnames)[0])[0]

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf existing res ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not path.exists(slf_file):
        raise TelemacException(\
               '... the provided slf_file does not seem to exist: '
               '{}\n\n'.format(slf_file))
    slf = Selafin(slf_file)
    slf.set_kd_tree()
    slf.set_mpl_tri()

    print('   +> support extraction')
    # Extract triangles and weigths in 2D
    support2d = []
    if showbar:
        ibar = 0
        pbar = ProgressBar(maxval=len(xys)).start()
    for xyi in xys:
        support2d.append(xys_locate_mesh(xyi, slf.ikle2, slf.meshx, slf.meshy,
                                         slf.tree, slf.neighbours))
        if showbar:
            ibar += 1
            pbar.update(ibar)
    if showbar:
        pbar.finish()
    # Extract support in 3D
    support3d = list(zip(support2d, len(xys)*[range(slf.nplan)]))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes BND header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bnd = Selafin('')
    bnd.fole = {}
    bnd.fole.update({'hook':open(bnd_file, 'wb')})
    bnd.fole.update({'name':bnd_file})
    bnd.fole.update({'endian':">"})     # big endian
    bnd.fole.update({'float':('f', 4)})  # single precision

    # Meta data and variable names
    bnd.title = ''
    bnd.nbv1 = len(varnames)
    # /!\ ELEVATION has to be the first variable
    # (for possible vertical re-interpolation within TELEMAC)

    bnd.varnames = []
    bnd.varunits = []
    for var, unit in zip(varnames, varunits):
        new_var = var + (16-len(var))*" "
        new_unit = unit + (16-len(unit))*" "
        bnd.varnames.append(new_var)
        bnd.varunits.append(new_unit)

    bnd.nvar = bnd.nbv1
    bnd.varindex = range(bnd.nvar)

    # Sizes and mesh connectivity
    bnd.nplan = slf.nplan
    # Number of nodes per boundary element  (ndp2 in 2D and ndp3 in 3D)
    bnd.ndp2 = 2
    bnd.ndp3 = 4
    bnd.npoin2 = len(bor)
    bnd.npoin3 = bnd.npoin2*slf.nplan
    bnd.iparam = [0, 0, 0, 0, 0, 0, bnd.nplan, 0, 0, 1]
    bnd.ipob2 = bor   # /!\ note that ipobo keeps the original numbering
    print('   +> masking and setting connectivity')
    # Set the array that only includes elements of geo.ikle2
    # with at least two nodes in bor
    array_1d = np.in1d(geo.ikle2, np.sort(bor-1))
    mask = geo.ikle2[np.where(np.sum(array_1d.reshape(geo.nelem2, geo.ndp2),
                                     axis=1) == 2)]
    # this ikle2 keeps the original numbering
    ikle2 = np.ravel(mask)[np.in1d(mask, np.sort(bor-1))].reshape(len(mask), 2)
    # ~~> re-numbering ikle2 as a local connectivity matrix
    knolg, _ = np.unique(np.ravel(ikle2), return_index=True)
    knogl = dict(zip(knolg, range(len(knolg))))
    bnd.ikle2 = - np.ones_like(ikle2, dtype=np.int)
    for k in range(len(ikle2)):
        # /!\ bnd.ikle2 has a local numbering, fit to the boundary elements
        bnd.ikle2[k] = [knogl[ikle2[k][0]], knogl[ikle2[k][1]]]
    # Last few numbers
    bnd.nelem2 = len(bnd.ikle2)
    if slf.nplan > 1:
        bnd.nelem3 = bnd.nelem2*(slf.nplan-1)
    else:
        bnd.nelem3 = bnd.nelem2
        bnd.ndp3 = bnd.ndp2
    # 3D structures
    if slf.nplan > 1:
        bnd.ipob3 = np.ravel(np.add(np.repeat(bnd.ipob2, slf.nplan)\
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
    bnd.datetime = slf.datetime
    bnd.tags['times'] = slf.tags['times']
    # VARIABLE extraction
    list_var = varnames[0]+": "
    for var in varnames[1:]:
        list_var += ";"+var+": "

    vrs = subset_variables_slf(list_var, slf.varnames)

    # Read / Write data, one time step at a time to support large files
    print('   +> reading / writing variables')
    if showbar:
        pbar = ProgressBar(maxval=len(slf.tags['times'])).start()
    zeros = np.zeros((bnd.npoin3, 1), dtype=np.float)
    for itime in range(len(slf.tags['times'])):
        data = get_value_history_slf(slf.file, slf.tags, [itime], support3d,
                                     slf.nvar, slf.npoin3, slf.nplan, vrs)
        data = np.reshape(np.transpose(np.reshape(np.ravel(data),
                                                  (bnd.nvar, bnd.npoin2,
                                                   bnd.nplan)),
                                       (0, 2, 1)),
                          (bnd.nvar, bnd.npoin3))
        bnd.append_core_time_slf(itime)
        bnd.append_core_vars_slf(data)
        if showbar:
            pbar.update(itime)
    if showbar:
        pbar.finish()

    # Close bnd_file
    bnd.fole['hook'].close()
