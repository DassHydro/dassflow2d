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
        Creates an atmospheric condition file from a global model (SLF form)
          by interpolating on a given GEO model domain.

    @history 16/12/2014 -- Sebastien E. Bourban
        Created from convertINI.py and createBND.py
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
from data_manip.conversion.convert_utm import to_latlon
from data_manip.formats.selafin import Selafin
from data_manip.extraction.parser_selafin import \
                  subset_variables_slf, get_value_history_slf
from pretel.meshes import xys_locate_mesh

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#
__author__ = "S.E.Bourban"
__date__ = "$16-Dec-2014 15:09:48$"

def generate_atm_parser(subparser):
    """
    Adding options for genreate_atm ot parser

    @param subparser (ArgumentParser) A parser

    @return (ArgumentParser) The updated parser
    """
    parser = subparser.add_parser('generate_atm',\
            help='Map weather type data (varying in space and time) contained'\
                 'into a SELAFIN, onto a SELAFIN file of your choosing'\
                 '(your MESH)')
    parser.add_argument(\
         "geo_file", default='',
         help="Geometry file on wich to map the weather type file")
    parser.add_argument(\
         "slf_file", default='',
         help="Weather data file from which to extract weather informations")
    parser.add_argument(\
         "atm_file", default='',
         help="Output file containing the weather information applied on the "
              "geometry file")
    parser.add_argument(\
            "--ll2utm", dest="ll2utm", default=None,
            help="assume outer file is in lat-long and open-bound file in UTM")

    return subparser

def generate_atm(geo_file, slf_file, atm_file, ll2utm):

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ slf new mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not path.exists(geo_file):
        raise TelemacException(\
            '... the provided geo_file does not '
            'seem to exist: {}\n\n'.format(geo_file))

# Find corresponding (x,y) in corresponding new mesh
    print('   +> getting hold of the GEO file')
    geo = Selafin(geo_file)
    if ll2utm is not None:
        zone = int(ll2utm[:-1])
        zone_letter = ll2utm[-1]
        x, y = to_latlon(geo.meshx, geo.meshy, zone, zone_letter)
    else:
        x = geo.meshx
        y = geo.meshy
    xys = np.vstack((x, y)).T

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
# ~~~~ writes ATM header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    atm = Selafin('')
    atm.fole = {}
    atm.fole.update({'hook':open(atm_file, 'wb')})
    atm.fole.update({'name':atm_file})
    atm.fole.update({'endian':">"})     # big endian
    atm.fole.update({'float':('f', 4)})  # single precision

    # Meta data and variable names
    atm.title = ''
    atm.varnames = []
    atm.varunits = []
    if 'WIND VELOCITY U ' in slf.varnames:
        atm.varnames.append('WIND VELOCITY U ')
        atm.varunits.append('M/S             ')
    if 'WIND VELOCITY V ' in slf.varnames:
        atm.varnames.append('WIND VELOCITY V ')
        atm.varunits.append('M/S             ')
    if 'SURFACE PRESSURE' in slf.varnames:
        atm.varnames.append('SURFACE PRESSURE')
        atm.varunits.append('UI              ')
    if 'AIR TEMPERATURE ' in slf.varnames:
        atm.varnames.append('AIR TEMPERATURE ')
        atm.varunits.append('DEGREES         ')
    if not atm.varnames:
        raise TelemacException(
            'There are no meteorological variables to convert!')
    atm.nbv1 = len(atm.varnames)
    atm.nvar = atm.nbv1
    atm.varindex = range(atm.nvar)

    # Sizes and mesh connectivity
    atm.nplan = slf.nplan          # it should be 2d but why the heack not ...
    atm.ndp2 = slf.ndp2
    atm.ndp3 = slf.ndp3
    atm.npoin2 = geo.npoin2
    atm.npoin3 = geo.npoin2*atm.nplan
    atm.nelem2 = geo.nelem2

    print('   +> setting connectivity')
    if atm.nplan > 1:
        atm.nelem3 = geo.nelem2*(atm.nplan-1)
        atm.ikle2 = geo.ikle2
        atm.ikle3 = \
            np.repeat(geo.npoin2*np.arange(atm.nplan-1),
                      geo.nelem2*atm.ndp3)\
              .reshape((geo.nelem2*(atm.nplan-1), atm.ndp3)) + \
            np.tile(np.add(np.tile(geo.ikle2, 2),
                           np.repeat(geo.npoin2*np.arange(2),
                                     geo.ndp2)),
                    (atm.nplan-1, 1))
        atm.ipob2 = geo.ipob2
        atm.ipob3 = np.ravel(np.add(np.repeat(geo.ipob2, atm.nplan)\
                                      .reshape((geo.npoin2, atm.nplan)),
                                    geo.npoin2*np.arange(atm.nplan)).T)
    else:
        atm.nelem3 = geo.nelem2
        atm.ikle2 = geo.ikle2
        atm.ikle3 = geo.ikle3
        atm.ipob2 = geo.ipob2
        atm.ipob3 = geo.ipob3
    atm.iparam = [0, 0, 0, 0, 0, 0, 0, np.count_nonzero(atm.ipob2), 0, 1]

    # Mesh coordinates
    atm.meshx = geo.meshx
    atm.meshy = geo.meshy

    print('   +> writing header')
    # Write header
    atm.datetime = slf.datetime
    atm.append_header_slf()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ writes ATM core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('   +> setting variables')
    # TIME and DATE extraction
    atm.tags['times'] = slf.tags['times']
    # VARIABLE extraction
    vrs = subset_variables_slf(';'.join([var+': ' for var in atm.varnames]),
                               slf.varnames)

    # Read / Write data, one time step at a time to support large files
    pbar = ProgressBar(maxval=len(slf.tags['times'])).start()
    for time in range(len(slf.tags['times'])):

        data = get_value_history_slf(slf.file, slf.tags, [time], support3d,
                                     slf.nvar, slf.npoin3, slf.nplan, vrs)
        # special cases ?
        atm.append_core_time_slf(time)
        atm.append_core_vars_slf(np.reshape(np.transpose(\
                  np.reshape(np.ravel(data),
                             (atm.nvar, atm.npoin2, atm.nplan)),
                  (0, 2, 1)),
                                            (atm.nvar, atm.npoin3)))
        pbar.update(time)
    pbar.finish()

    # Close atm_file
    atm.fole['hook'].close()
