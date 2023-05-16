r"""@author TELEMAC-MASCARET Consortium

    @brief
    Tools for handling boundary conditions files
"""

from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import numpy as np

from utils.files import get_file_content, put_file_content


class Conlim(object):
    """ Boundary conditions class """

    def __init__(self, file_name):
        """
        @brief Initializes the Conlim object

        @param file_name (string): name of the boundary conditions file

        @return None
        """

        self.file_name = file_name
        # ~~> Columns of integers and floats
        dtype = np.dtype = [('lih', '<i4'), ('liu', '<i4'), ('liv', '<i4'),
                            ('h', '<f4'), ('u', '<f4'), ('v', '<f4'),
                            ('au', '<f4'), ('lit', '<i4'), ('t', '<f4'),
                            ('at', '<f4'), ('bt', '<f4'), ('n', '<i4'),
                            ('line', '<i4')]

        if file_name != '':
            # ~~> Number of boundary points ( tuple() necessary for dtype
            # parsing )
            core = [tuple(line.strip().split()[0:13])
                    for line in get_file_content(file_name)]
            self.nptfr = len(core)
            self.bor = np.array(core, dtype)
            # ~~> Dictionary of kfrgl
            self.kfrgl = dict(zip(self.bor['n']-1, range(self.nptfr)))

            # ~~> Filtering indices
            self.index = np.array(range(self.nptfr), dtype=np.int)
        else:
            self.nptfr = 0
            self.bor = None
            self.kfrgl = None
            self.index = None
        self.nptir = 0
        self.ifapar = None
        self.por = None
        self.nfrliq = 0

    def set_numliq(self, closed_contours):
        """
        @brief Sets the number of liquid boundaries

        @param closed_contours (list): list of closed contours

        @return None
        """

        # ~~> Part linkages
        self.nptir = {}
        self.ifapar = {}

        # ~~> Extra columns
        dtype = np.dtype = [('is', '<i4'), ('xs', '<f4'),
                            ('ys', '<f4'), ('lq', '<i4')]
        core = [tuple(' 0 0.0 0.0 0'.split()) for _ in range(self.nptfr)]
        self.por = np.array(core, dtype)

        # ~~> Initiates NUMLIQ
        solids = self.bor['lih'] == 2
        self.nfrliq = 0  # liquid boundary numbers start at 1

        # ~~> Finds TELEMAC's south-east corner

        # ~~> Counts NUMLIQ
        # for the domain boundary and eavery islands
        for contour in closed_contours:
            # ~~> look for a solid point
            inode = 0
            while not solids[self.kfrgl[contour[inode]]] and \
                    inode < len(contour):
                inode += 1
            if inode == len(contour):  # No solid boundary found in contour
                self.nfrliq += 1
                for i in contour:
                    self.por['lq'][self.kfrgl[i]] = self.nfrliq
                continue  # go to next contour
            # ~~> list all liquid boundaries on contour
            wassolid = True
            for i in contour:  # this makes sure to go around contour only once
                isolid = self.kfrgl[contour[inode % len(contour)]]
                if not solids[isolid]:
                    if wassolid:
                        self.nfrliq += 1
                        wassolid = False
                    self.por['lq'][isolid] = self.nfrliq
                else:
                    wassolid = True
                inode += 1

        return

    def put_content(self, file_name):
        """
        @brief Writes the liquid boundaries file

        @param file_name (string): name of the liquid boundaries file
        """

        # standard part of the Conlim
        core = []
        for ifr in range(self.nptfr):
            if self.index[ifr] != 0:
                line = (' '.join(['{0[' + repr(i) + ']}'
                        for i in range(len(self.bor[ifr]))]))\
                     .format(self.bor[ifr])
                if self.nptir != {} and self.ifapar != {}:
                    line += ' '+repr(self.por[ifr][0])
                    line += ' '+repr(float(self.por[ifr][1]))
                    line += ' '+repr(float(self.por[ifr][2]))
                    line += ' '+repr(self.por[ifr][3])
                core.append(line)

        if self.nptir != {} and self.ifapar != {}:
            # parrallel part 1 of the Conlim -- format with i7 /!\
            ntr = self.nptir.keys()
            ntr.sort()
            core.append(str(len(ntr)))
            for itr in ntr:
                line = self.nptir[itr]
                line.sort()
                line.extend([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
                s_line = repr(itr+1).rjust(7) + ' ' + \
                    repr(line[0]).rjust(7) + ' ' + \
                    repr(line[1]).rjust(7) + ' ' + \
                    repr(line[2]).rjust(7) + ' ' + \
                    repr(line[3]).rjust(7) + ' ' + \
                    repr(line[4]).rjust(7) + ' ' + \
                    repr(line[5]).rjust(7) + ' ' + \
                    repr(line[6]).rjust(7) + ' ' + \
                    repr(line[7]).rjust(7) + ' ' + \
                    repr(line[8]).rjust(7)
                core.append(s_line)
        # parrallel part 2 of the Conlim -- format with i9 /!\
            ntr = self.ifapar.keys()
            ntr.sort()
            core.append(str(len(ntr)))
            for itr in ntr:
                line = self.ifapar[itr]
                s_line = repr(itr+1).rjust(9) + ' ' + \
                    repr(line[1]).rjust(9) + ' ' + \
                    repr(line[3]).rjust(9) + ' ' + \
                    repr(line[5]).rjust(9) + ' ' + \
                    repr(line[0]+1).rjust(9) + ' ' + \
                    repr(line[2]+1).rjust(9) + ' ' + \
                    repr(line[4]+1).rjust(9)
                core.append(s_line)

        core.append("")
        put_file_content(file_name, core)
