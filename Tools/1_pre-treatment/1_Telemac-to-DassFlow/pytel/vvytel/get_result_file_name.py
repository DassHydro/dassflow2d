r"""@author TELEMAC-MASCARET Consortium

    @note ... this work is based on a collaborative effort between
                 TELEMAC-MASCARET consortium members
    @brief Function for validation of the Python API
"""
from os import path, environ
from execution.telemac_cas import TelemacCas
###############################################################################
KEY_RES = {'telemac2d':'RESULTS FILE',
           'telemac3d':'3D RESULT FILE',
           'tomawac':'2D RESULTS FILE',
           'artemis':'RESULTS FILE'}
###############################################################################
def get_result_file_name(module, cas_file):
    """
       Returns the name of the result file for a given case.

       @param module	name of the telemac module
       @param cas_file	name of the telemac steering file
    """
    dico_file = path.join(environ['HOMETEL'], 'sources', module, module+'.dico')
    cas = TelemacCas(cas_file, dico_file)
    res_file = cas.get(KEY_RES[module])
    return res_file

