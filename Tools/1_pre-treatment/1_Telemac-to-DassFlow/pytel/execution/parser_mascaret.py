r"""@author Christophe Coulet

    @history 27/06/2017 -- Christophe Coulet:
            Creation of parser of XCAS file
            Return the list of Mascaret Input File

    @history 20/07/2017 -- Christophe Coulet:
            Rename the parser which become parserMascaret
            Adding capabilities to read Opthyca and Rubens results


    @brief
"""
# _____             ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import xml.etree.ElementTree as ET
from os import path
# ~~> dependencies towards other pytel/modules
from config import CFGS
from execution.telemac_cas import TelemacCas


def scan_xcas(fle):
    """
    @brief : read the xml file to extract the list of input file
    :param fle: xcas file of mascaret computation
    :return: list of file needed for computation
    """
    inputfile = []
    tree = ET.parse(fle)
    root = tree.getroot()
    root2 = root[0]

    # looking for geometry
    inputfile.append(root2.find('parametresGeometrieReseau')
                     .find('geometrie')
                     .find('fichier').text)

    # looking for laws
    lois = root2.find('parametresLoisHydrauliques').find('lois')
    for loi in lois:
        inputfile.append(loi.find('donnees').find('fichier').text)

    # looking for initial conditions
    linits = root2.find('parametresConditionsInitiales').find('ligneEau')
    if linits.find('LigEauInit').text == 'true':
        inputfile.append(linits.find('fichLigEau').text)

    # looking for "casier"
    if root2.find('parametresCasier') is not None:
        inputfile.append((root2.find('parametresCasier')
                               .find('fichierGeomCasiers').text))

    # looking for "paramtresPhysique"
    if root2.find('parametresTraceur') is not None:
        root_tracer = root2.find('parametresTraceur')
        if root_tracer.find(
                    'parametresConcentrationsInitialesTraceur') is not None:
            inputfile.append(
                    root_tracer.find(
                        'parametresConcentrationsInitialesTraceur')
                    .find('fichConcInit').text)

        if root_tracer.find('parametresNumeriquesQualiteEau') is not None:
            inputfile.append(root_tracer.find(
                'parametresNumeriquesQualiteEau')
                                        .find('fichParamPhysiqueTracer').text)
            inputfile.append(root_tracer.find(
                'parametresNumeriquesQualiteEau')
                                        .find('fichMeteoTracer').text)

        lois = root_tracer.find('parametresLoisTraceur').find('loisTracer')
        for loi in lois:
            inputfile.append(loi.find('fichier').text)

    # looking for "Courlis"
    if root2.find('parametresGeneraux').find('optionCourlis') is not None:
        inputfile.append((root2.find('parametresGeneraux')
                         .find('fichierMotCleCourlis').text))

        casfile = root2.find('parametresGeneraux')\
                       .find('fichierMotCleCourlis').text
        print(casfile)
        dicofile = path.join(CFGS.get_root(), "sources", "mascaret", "data",
                             "dico_Courlis.txt")

        cas = TelemacCas(casfile, dicofile)
        geo_courlis = cas.get('FICHIER DE GEOMETRIE COURLIS')
        inputfile.append(geo_courlis)

    return inputfile
