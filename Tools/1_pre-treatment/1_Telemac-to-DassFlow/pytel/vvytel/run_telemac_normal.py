r"""@author TELEMAC-MASCARET Consortium

    @note ... this work is based on a collaborative effort between
                 TELEMAC-MASCARET consortium members
    @brief
"""

from os import system
###############################################################################
def run_telemac_normal(module, cas, ncsize):
    """
    Normal run of Telemac
    @param module	Name of the module
	@param cas		The name of the steering file
    @param ncsize	Number of parallel processors
    """
    cmd = module+'.py --ncsize='+str(ncsize)+' '+cas+' > run.log'
    print(cmd)
    system(cmd)
    passed = False
    with open('run.log', 'r') as fobj:
        for line in fobj.readlines():
            if "My work is done" in line:
                passed = True
    return passed

