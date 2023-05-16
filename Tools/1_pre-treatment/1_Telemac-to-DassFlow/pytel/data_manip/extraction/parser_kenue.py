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
            Tools for handling Kenue native files in python.
            Kenue and its related software (Blue Kenue, Greem Kenue, )
            are property of the NRC Canadian Hydrualics Centre

    @details
            Contains getI2S, getI3S and putI2S, putI3S, which read/write
            python variables into ASCII I2S and I3S files

    @history 26/12/2011 -- Sebastien E. Bourban:
            First trial at parsing I2S and I3S

    @history 09/01/2012 -- Sebastien E. Bourban:
            Addition of XY and XYZ parsing

    @history 13/01/2012 -- Sebastien E. Bourban:
            Creates InS class with associated methods including:
            + removeDuplicates (remove duplicated points based on proximity)
            + make_clockwise (make closed loops clockwise)
            + make_anti_clockwise (make closed loops anti-clockwise)
            + smoothSubdivise (add points and weigthed average move)
            + smoothSubsampleDistance (remove points based on proximity)
            + smoothSubsampleAngle (remove points based on flatness)

    @history 29/02/2012 -- Sebastien E. Bourban:
            The type of a contour can now be checked automatically so that
            mixed set can be stored into one I2S/I3S file. To use this
            feature, type must be absent from the putInS call.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import re
import numpy as np
from os import path
# ~~> dependencies towards other pytel/modules
from utils.files import get_file_content, put_file_content
from utils.geometry import is_close
from utils.exceptions import TelemacException

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#
KEN_HEADER = re.compile(r'[#:]')

ASC_FILE_TYPE = re.compile(r':FileType\s(?P<type>\b\w\w\w\b)') #\s(?P<after>.*?)\Z')
ASC_ATTRIBUTE_NAME = re.compile(r':AttributeName\s(?P<number>\b(|[^a-zA-Z(,])(?:(\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')

VAR_1INT = re.compile(r'(?P<before>[^\'"]*?)\b(?P<number>[+-]?(|[^a-zA-Z(,])(?:(\d+)(\b|[^a-zA-Z,)])))(?P<after>[^.].*?)\Z')
VAR_1DBL = re.compile(r'(?P<number>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<after>.*?)\Z')
VAR_2DBL = re.compile(r'(?P<number1>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<number2>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')
VAR_3DBL = re.compile(r'(?P<number1>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<number2>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<number3>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')
VAR_1STR = re.compile(r'(?P<string>)(".*?")[\s,;]*(?P<after>.*?)\Z')

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

# _____                              _______________________________
# ____/ Principal Class for I2S/I3S /______________________________/
#


class InS(object):

    def __init__(self, file_name):
        self.atrbut = None
        self.natrbut = 0
        self.oatrbut = None
        self.file_type = 0
        self.head = None
        self.poly = None
        self.vals = None
        self.type = None
        self.npoin = None
        self.npoly = 0

        # file parsing is based on the name of the extension
        _, tail = path.splitext(file_name)
        # ~~> Case of a Kenue type i2s/i3s file
        if tail in ['.i2s', '.i3s']:
            self.parser_content(file_name)
        else:
            raise TelemacException(
                    '\nThe polygon file extension is required to be '
                    'either i2s or i3s')

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Parse Content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    def parser_content(self, file_name):
        # TODO: Read the whole header

        # ~~ Get all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        core = get_file_content(file_name)

        # ~~ Parse head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        icore = 0
        self.atrbut = {}
        self.natrbut = 0
        self.oatrbut = []
        while re.match(KEN_HEADER, core[icore].strip()):
            # ~~> instruction FileType
            proc = re.match(ASC_FILE_TYPE, core[icore].strip())
            if proc:
                self.file_type = proc.group('type').lower()
            # ~~> instruction AttributeName
            proc = re.match(ASC_ATTRIBUTE_NAME, core[icore].strip())
            if proc:
                self.natrbut += 1
                if self.natrbut == int(proc.group('number')):
                    self.oatrbut.append(proc.group('after').strip())
                    self.atrbut.update({self.oatrbut[-1]: []})
                else:
                    raise TelemacException(
                            '... Could not read the order of your Attributes: '
                            '{}'.format(core[icore]))
            # ... more instruction coming ...
            icore += 1
        self.head = core[0:icore]
        # /!\ icore starts the body

        # ~~ Parse body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # This is also fairly fast, so you might not need a progress bar
        self.poly = []
        self.vals = []
        self.type = []
        self.npoin = 0
        while icore < len(core):
            if core[icore].strip() == '':
                icore += 1
                continue
            # ~~> polygon head
            proc = re.match(VAR_1INT, core[icore].strip())
            if not proc:
                raise TelemacException(
                        '\nCould not parse the following polyline header: '
                        '{}'.format(core[icore].strip()))
            nrec = int(proc.group('number'))
            after = proc.group('after').strip().split()
            if len(after) != self.natrbut:
                if self.natrbut != 0:
                    raise TelemacException(
                        '... Could not find the correct number of attribute:'
                        '{}, {} expected'.format(core[icore].strip(),
                                                 self.natrbut))
                else:
                    self.natrbut = len(after)
                    self.oatrbut = range(1, len(after)+1)
                    self.atrbut = \
                        dict([(i+1, [after[i]]) for i in range(len(after))])
            else:
                for i in range(len(self.oatrbut)):
                    self.atrbut[self.oatrbut[i]].append(after[i])
            xyi = []
            val = []
            icore += 1
            for irec in range(nrec):
                nbres = core[icore+irec].strip().split()
                proc = re.match(VAR_1DBL, nbres[0])
                if not proc:
                    proc = re.match(VAR_1INT, nbres[0])
                    if not proc:
                        raise TelemacException(
                            '\nCould not parse the following polyline record: '
                            '{}'.format(core[icore+irec].strip()))
                nbres[0] = float(proc.group('number'))
                procd = re.match(VAR_1DBL, nbres[1])
                proci = re.match(VAR_1INT, nbres[1])
                if procd:
                    nbres[1] = float(procd.group('number'))
                elif proci:
                    nbres[1] = float(procd.group('number'))
                xyi.append(nbres[:2])
                val.append(nbres[2:])
            if xyi != []:
                cls = 0
                if is_close(xyi[0], xyi[len(xyi)-1], size=10):
                    xyi.pop(len(xyi)-1)
                    val.pop(len(val)-1)
                    cls = 1
                self.poly.append(np.asarray(xyi, dtype=np.float))
                self.vals.append(np.asarray(val, dtype=np.float))
                self.type.append(cls)
            self.npoin += len(xyi)
            icore += nrec

        self.npoly = len(self.poly)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Write-up Content ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Note:
    #    + file parsing is based on the name of the extension
    #
    def put_content(self, file_name, head=None):

        # ~~> file extension processing
        _, tail = path.splitext(file_name)
        if tail[1:] != self.file_type:
            if head is not None:
                head = ['\n'.join(head).replace(':FileType '+self.file_type,
                                                ':FileType '+tail[1:])]
            self.file_type = tail[1:]

        # ~~> write head
        if head != []:
            core = head
        else:
            core = [':FileType '+self.file_type+' ASCII EnSim 1.0',
                    ':Application BlueKenue',
                    ':Version 3.2.24',
                    ':WrittenBy sebourban',
                    ':CreationDate Thu, Dec 08, 2011 02:47 PM',
                    ':Name ' + path.basename(file_name),
                    # ':AttributeName 1 level',
                    # ':AttributeType 1 float',
                    # ':AttributeUnits 1 m',
                    ':EndHeader']

        # ~~> look for closed lines
        if self.type == []:
            for i_p in self.poly:
                if is_close(i_p[0][:2], i_p[len(i_p)-1][:2]):
                    self.type.append(1)
                else:
                    self.type.append(0)

        # ~~> fill-up empty attributes
        if self.atrbut == {}:
            self.atrbut = {1: ['ArbitraryName1']}
            for _ in self.poly:
                self.atrbut[1].append(0)
            self.oatrbut = [1]

        # ~~> fill-up attribute names
        if head == []:
            for i, val in enumerate(self.oatrbut):
                core.insert(-1, ':AttributeName '+repr(i+1)+' '+repr(val))

        # ~~ Write body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i, ipoly, ival, itype in zip(range(len(self.poly)),
                                         self.poly,
                                         self.vals,
                                         self.type):
            i_l = len(ipoly)
            if i_l != 0 and not is_close(ipoly[0], ipoly[len(ipoly)-1]):
                i_l += itype
            line = repr(i_l)
            for oatr in self.oatrbut:
                line += ' ' + str(self.atrbut[oatr][i])
            core.append(line)
            if self.file_type == 'i2s':
                for xyi in ipoly:
                    core.append(repr(xyi[0])+' '+repr(xyi[1]))
                if i_l != len(ipoly):
                    core.append(repr(ipoly[0][0])+' '+repr(ipoly[0][1]))
            elif self.file_type == 'i3s':
                if np.shape(ival)[1] == 0:
                    for xyi in ipoly:
                        core.append(repr(xyi[0])+' '+repr(xyi[1])+' 0.0')
                    if i_l != len(ipoly):
                        core.append(repr(ipoly[0][0])+' ' +
                                    repr(ipoly[0][1])+' 0.0')
                else:
                    for xyi, val in zip(ipoly, ival):
                        core.append(repr(xyi[0])+' '+repr(xyi[1])+' ' +
                                    ' '.join([repr(v) for v in val]))
                    if i_l != len(ipoly):
                        core.append(repr(ipoly[0][0])+' '+repr(ipoly[0][1]) +
                                    ' '+' '.join([repr(v) for v in ival[0]]))

        # ~~ Put all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        put_file_content(file_name, core)

# _____                  ___________________________________________
# ____/ Toolbox for XYZ /__________________________________________/
#


def get_xyn(file_name):
    # TODO: Read the whole header, for the time being head is copied
    #       over
    # TODO: Read multiple variables depending on type and on a list

    # ~~ Get all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core = get_file_content(file_name)

    # ~~ Parse head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    icore = 0
    file_type = None
    while re.match(KEN_HEADER, core[icore]):
        # ~~> instruction FileType
        proc = re.match(ASC_FILE_TYPE, core[icore])
        if proc:
            file_type = proc.group('type').lower()
        # ... more instruction coming ...
        icore += 1
    head = core[0:icore]
    if file_type is None:
        proc = re.match(VAR_3DBL, core[icore]+' ')
        if not proc:
            proc = re.match(VAR_2DBL, core[icore]+' ')
            if not proc:
                raise TelemacException(
                        '\nCould not parse the first record: '
                        '{}'.format(core[icore]))
            else:
                file_type = 'xy'
        else:
            file_type = 'xyz'

    # /!\ icore starts the body

    # ~~ Parse body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is also fairly fast, so you might not need a progress bar
    xyz = []  # ; pbar = ProgressBar(maxval=len(core)).start()
    while icore < len(core):
        if file_type == 'xy':
            proc = re.match(VAR_2DBL, core[icore]+' ')
            if not proc:
                raise TelemacException(
                        '\nCould not parse the following xyz record: '
                        '{}'.format(core[icore]))
            xyz.append([float(proc.group('number1')),
                        float(proc.group('number2'))])
        elif file_type == 'xyz':
            proc = re.match(VAR_3DBL, core[icore]+' ')
            if not proc:
                raise TelemacException(
                        '\nCould not parse the following xyz record: '
                        '{}'.format(core[icore]))
            xyz.append([float(proc.group('number1')),
                        float(proc.group('number2')),
                        float(proc.group('number3'))])
        icore += 1
    # pbar.finish()

    return head, file_type, xyz


def put_xyn(fle, head, file_type, xyz):

    # ~~ Write head ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core = head
    # if head != []: core = head
    # <else: core = [':FileType '+file_type+' ASCII EnSim 1.0',
    #   ':Application BlueKenue', ':Version 3.2.24',
    #   ':WrittenBy sebourban', ':CreationDate Thu, Dec 08, 2011 02:47 PM',
    #   ':Name ' + path.basename(file_name),
    #   ':AttributeUnits 1 m',
    #   ':EndHeader' ]

    # ~~ Write body ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if file_type == 'xy':
        for _ in xyz:
            core.append(str(xyz[0])+' '+str(xyz[1]))
    elif file_type == 'xyz':
        for _ in xyz:
            core.append(str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2]))

    # ~~ Put all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    put_file_content(fle, core)

    return
