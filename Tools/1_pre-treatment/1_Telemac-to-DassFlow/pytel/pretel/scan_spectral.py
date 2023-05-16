r"""@author Sebastien E. Bourban

"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
from pretel.scan_selafin import ScanSelafin
import numpy as np

class ScanSpectral(ScanSelafin):
    """
    Dump of spectrum file information
    """

    def print_header(self):
        """
        Printing header info
        """
        ScanSelafin.print_header(self)
        ndir = self.slf.npoin3-self.slf.nelem3
        nfreq = self.slf.npoin3//ndir
        print("FREQUENCIES "+str(nfreq)+" / min: [ "+str(self.slf.meshy[0])+\
                "]  / max: [ "+str(self.slf.meshy[self.slf.npoin3-ndir])+"]")
        print(' '.join([repr(i) for i in self.slf.meshy[0::ndir]]))
        print("DIRECTIONS  "+str(ndir)+" / angle: "+str(360./ndir))

    def print_core(self, na):
        """
        Print variable information
        """
        npoin3 = self.slf.npoin3
        nelem3 = self.slf.nelem3
        for v in range(self.slf.nbv1):
            print("VARIABLE     : "+self.slf.varnames[v])
            for itime in range(len(self.slf.tags['times'])):
                varsor = self.slf.get_variables_at(itime,
                                                   [self.slf.varindex[v]])
                print("    / TIME: "+str(self.slf.tags['times'][itime]))
                # na significant figures
                accuracy = np.power(10.0,
                                    -na+np.floor(np.log10(abs(np.max(varsor)))))
                print('\nFACTOR '+ str(accuracy))
                spes = np.reshape(np.array((varsor/accuracy), dtype=np.int),
                                  (npoin3//(npoin3-nelem3), -1)) # by frequency
                for spe in spes:
                    print(' '.join([repr(i).rjust(na+1) for i in spe]))
                print('\n')
        for v in range(self.slf.nbv2):
            print("CLANDESTINE  : "+str(self.slf.cldnames[v]))
            for itime in range(len(self.slf.tags['times'])):
                varsor = self.slf.get_variables_at(\
                          itime, [self.slf.varindex[v+self.slf.nbv1]])
                print("    / TIME: "+str(self.slf.tags['times'][itime]))
                # na significant figures
                accuracy = np.power(10.0,
                                    -na+np.floor(np.log10(abs(np.max(varsor)))))
                print('\nFACTOR '+ str(accuracy))
                spes = np.reshape(np.array((varsor/accuracy), dtype=np.int),
                                  (npoin3//(npoin3-nelem3), -1)) # by frequency
                for spe in spes:
                    print(' '.join([repr(i).rjust(na+1) for i in spe]))
                print('\n')
