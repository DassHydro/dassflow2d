r"""@author Telemac-Mascaret Consortium

    @brief
            Tools for handling Time Series files (such as LQD) in python.

    @details
            Contains get_lqd and put_lqd, which read/write python variables
            into ASCII LQD files

    @history 26/12/2011 -- Sebastien E. Bourban:
            First trial at writting LDQ having sampled SELAFIN files at nodes
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import re
import numpy as np
# ~~> dependencies towards other pytel/modules
from utils.files import get_file_content, put_file_content
from utils.exceptions import TelemacException

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#
LQD_HEADER = re.compile(r'#')

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

# _____                  ___________________________________________
# ____/ Toolbox for LDQ /__________________________________________/
#


def get_lqd(file_name):
    # ~~ Get all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core = get_file_content(file_name)

    # ~~ Parse head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    icore = 0
    while re.match(LQD_HEADER, core[icore]):
        icore += 1
    head = core[0:icore]
    # /!\ icore starts the body

    # ~~ Parse variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    jcore = icore+1
    while icore < len(core) and jcore < len(core):
        if re.match(LQD_HEADER, core[icore]):
            icore += 1
            jcore += 1
            continue
        if re.match(LQD_HEADER, core[jcore]):
            jcore += 1
            continue
        core[icore].replace(',', ' ')
        core[jcore].replace(',', ' ')
        # ~~> Variable header
        if core[icore].split()[0].upper() != 'T':
            raise TelemacException(
                    '\nThought I would find T for this LQD file on '
                    'this line: {}'.format(core[icore]))
        if len(core[icore].split()) != len(core[jcore].split()):
            raise TelemacException(
                '\nThought I to corresponding units for this LQD file on '
                'this line: {}'.format(core[jcore]))
        vrs = zip(core[icore].upper().split(), core[jcore].upper().split())

    # ~~ Size valid values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    icore = jcore+1
    while icore < len(core):
        if not re.match(LQD_HEADER, core[jcore]):
            icore += 1

    # ~~ Parse body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is also fairly fast, so you might not need a progress bar
    time = np.zeros(jcore-icore)
    z = np.zeros(len(vrs)-1, jcore-icore)
    itime = 0
    for icore in core[jcore+1:]:
        if re.match(LQD_HEADER, icore):
            continue
        values = icore.replace(',', ' ').split()
        time[itime] = float(values[0])
        for ivar in range(len(values[1:])):
            z[itime][ivar] = float(values[ivar])

    return head, vrs[1:], time, z


def put_lqd(fle, head, vrs, date0, time, xyz):

    # ~~ Write head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core = head
    from datetime import datetime
    core.append("#Date&Time: "+datetime(*date0).strftime("%d/%m/%Y %H:%M:%S"))

    names = 'T'
    units = 's'
    if xyz.ndim == 2:
        for name, unit in vrs:
            names += ' ' + name.strip().replace(' ', '_')
            units += ' ' + unit.strip().replace(' ', '_')
        core.append(names+'\n'+units)
    elif xyz.ndim == 3:
        for ivar in range(len(vrs[0])):
            for inod in vrs[1]:
                names += ' ' + vrs[0][ivar][0].strip().replace(' ', '_') + \
                            '(' + str(inod) + ')'
                units += ' ' + vrs[0][ivar][1].strip().replace(' ', '_')
        core.append(names+'\n'+units)
    if xyz.ndim == 2:
        for itim in range(xyz.shape[1]):
            line = str(time[itim])
            for ivar in range(xyz.shape[0]):
                line += ' ' + str(xyz[ivar][itim])
            core.append(line)
    elif xyz.ndim == 3:
        for itim in range(xyz.shape[2]):
            line = str(time[itim])
            for ivar in range(xyz.shape[0]):
                for inod in range(xyz.shape[1]):
                    line += ' ' + str(xyz[ivar][inod][itim])
            core.append(line)

    # ~~ Put all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    put_file_content(fle, core)

    return


class LQD(object):
    file_name = ''
    head = []

    def __init__(self, file_name='', vrs=None, date=None,
                 times=None, series=None):
        if file_name != '':  # read from file
            self.head, self.vrs, self.times, self.series = \
                    get_lqd(self.file_name)
        else:              # set content values
            self.vrs = vrs
            self.times = times
            self.series = series
        if date is None:
            self.date = [1972, 7, 13, 17, 24, 27]
        else:
            self.date = date

    def put_content(self, file_name):
        put_lqd(file_name, self.head, self.vrs, self.date,
                self.times, self.series)
