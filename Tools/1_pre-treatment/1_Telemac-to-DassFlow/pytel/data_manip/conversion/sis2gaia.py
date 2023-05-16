#!/usr/bin/python3
"""
Function to convert a sisyphe steering file into a gaia steering file
"""

from os import path, environ, remove
from execution.telemac_cas import TelemacCas


def sis2gaia_parser(subparser):
    """
    Defines the arguments for argparse

    @param subparser (ArgumentParser) The parser

    @returns (ArgumentParser)  The updated parser
    """
    parser = subparser.add_parser("sis2gaia",
                                  help="Conversion of the steering\
                                   file from SISYPHE to GAIA")
    parser.add_argument('sis_cas', help='Sisyphe steering file')
    parser.add_argument('gaia_cas', help='Gaia steering file')

    return subparser


def sis2gaia(sis_cas_file, gaia_cas_file):
    """
    Convert a sisyphe steering file into a gaia steering file

    @param sis_cas_file (string) Sisyphe cas file
    @param gaia_cas_file (string) Gaia cas file
    """
    trans_key = {
      'NUMBER OF BED MODEL LAYERS':
      'NUMBER OF LAYERS FOR INITIAL STRATIFICATION',
      'SOLVER FOR SUSPENSION': 'SOLVER FOR DIFFUSION OF SUSPENSION',
      'SOLVER OPTION FOR SUSPENSION':
      'SOLVER OPTION FOR DIFFUSION OF SUSPENSION',
      'PRECONDITIONING FOR SUSPENSION':
      'PRECONDITIONING FOR DIFFUSION OF SUSPENSION',
      'SOLVER ACCURACY FOR SUSPENSION': 'ACCURACY FOR DIFFUSION OF SUSPENSION',
      'SUSPENSION': 'SUSPENSION FOR ALL SANDS',
      'REFERENCE CONCENTRATION FORMULA':
      'SUSPENSION TRANSPORT FORMULA FOR ALL SANDS',
      'TETA SUSPENSION': 'TETA IMPLICITATION FOR SUSPENSION',
      'CRITICAL SHEAR VELOCITY FOR MUD DEPOSION':
      'CLASSES CRITICAL SHEAR STRESS FOR MUD DEPOSITION',
      'BED LOAD': 'BED LOAD FOR ALL SANDS',
      'BED-LOAD TRANSPORT FORMULA': 'BED-LOAD TRANSPORT FORMULA FOR ALL SANDS',
      'SEDIMENT DIAMETERS': 'CLASSES SEDIMENT DIAMETERS',
      'SETTLING VELOCITIES': 'CLASSES SETTLING VELOCITIES',
      'SEDIMENT DENSITY': 'CLASSES SEDIMENT DENSITY',
      'HIDING FACTOR FOR PATICULAR SIZE CLASS': 'CLASSES HIDING FACTOR',
      'SHIELDS PARAMETERS': 'CLASSES SHIELD PARAMETERS',
      'INITIAL FRACTION FOR PARTICULAR SIZE CLASS':
      'CLASSES INTITIAL FRACTION',
      'PARTHENIADES CONSTANT': 'LAYERS PARTHENIADES CONSTANT',
      'MUD CONCENTRATION PER LAYER': 'LAYERS MUD CONCENTRATION',
      'CRITICAL EROSION SHEAR STRESS OF THE MUD':
      'LAYERS CRITICAL EROSION SHEAR STRESS OF THE MUD',
      'MASS TRANSFER PER LAYER': 'LAYERS MASS TRANSFER',
      'NON COHESIVE BED POROSITY': 'LAYERS NON COHESIVE BED POROSITY'
                }
    removed_keys = [
      'STATIONARY MODE',
      'CONSTANT FLOW DISCHARGE',
      'NUMBER OF ITERATIONS FOR TELEMAC',
      'CRITERION TO UPDATE THE FLOW',
      'CRITICAL EVOLUTION RATIO',
      'NUMBER OF TIME STEP',
      'TIME STEP',
      'OPTION FOR THE TREATMENT OF NON ERODABLE BEDS',
      'GRAIN FEEDING',
      'MASS CONCENTRATION',
      'MIXED SEDIMENTS',
      'MEAN DIAMETER OF THE SEDIMENT',
      'HYDRODYNAMIC FILE',
      'D90',
      'FORMALATION FOR DEPOSITION AND EROSION',
      'STARTING TIME IF THE HYDROGRAM',
      'NUMBER OF TIDES OR FLOODS',
      'TIDE PERIOD',
      'PRECONDITIONING',
      'SOLVER',
      'MAXIMUM OF ITERATIONS FOR SOLVER',
      'SOLVER OPTION',
      'SOLVER ACCURACY',
      'PARTITIONNING TOOL',
      'WATER DENSITY',
      'GRAVITY ACCELERATION',
      'FRICTION COEFFICIENT',
      'LAW OF BOTTOM FRICTION',
      'DIFFUSION',
      'NUMBER IF SIZE-CLASSES OF BED MATERIAL',
      'COHESIVE SEDIMENTS',
      'VERTICAL GRAIN SORTING MODEL',
      'C-VSM MAIXMUM SECTIONS',
      'C-VSM FULL PRINTOUT PERIOD',
      'C-VSM PRINTOUT SELECTION',
      'C-VSM DYNAMIC ALT MODEL',
      'MASS-LUMPING',
      'TYPE OF ADVECTION',
      'SUPG OPTION',
      'OPTION FOR THE DIFFUSION OF TRACER',
      'NUMBER OF CORRECTIONS OF DISTRIBUTIVE SCHEMES',
      'NUMBER OF SUB-STEPS OF DISTRIBUTOVE SCHEMES',
      'TREATMENT OF FLUXES AT THE BOUNDARIES',
      'INITIAL SUSPENSION CONCENTRATIONS',
      'OPTION FOR THE DISPERSION',
      'SCHEME OPTION FOR ADVECTION',
      'CONSOLIDATION MODEL',
      'GEL CONCENTRATION',
      'MAXIMUM CONCENTRATION',
      'PERMEABILITY COEFFICIENT',
      'MUD CONSOLIDATION',
      'TETA'
                   ]

    sis_dico = path.join(environ['HOMETEL'], 'sources',
                         'sisyphe', 'sisyphe.dico')
    gaia_dico = path.join(environ['HOMETEL'], 'sources', 'gaia', 'gaia.dico')

    print(sis_cas_file, gaia_cas_file)
    sis_cas = TelemacCas(sis_cas_file, sis_dico)

    if path.exists(gaia_cas_file):
        remove(gaia_cas_file)

    gaia_cas = TelemacCas(gaia_cas_file, gaia_dico, access='w')

    for key, val in sis_cas.values.items():
        if key in removed_keys:
            # Special treatment
            pass
        else:
            # Simple translation
            new_key = trans_key.get(key, key)
            gaia_cas.set(new_key, val)

    gaia_cas.write(gaia_cas_file)
