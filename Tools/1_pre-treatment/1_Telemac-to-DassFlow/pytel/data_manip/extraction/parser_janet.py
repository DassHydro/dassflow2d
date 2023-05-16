r"""@author Sebastien E. Bourban and Michael S. Turnbull

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
            Tools for handling Janet native files in python.
            Janet and its related software (..., )
            are property of Smile Consulting

    @details
            Contains

    @history 26/12/2011 -- Sebastien E. Bourban:
            First trial at parsing I2S and I3S

    @history 13/01/2012 -- Sebastien E. Bourban:
            Creates Insel class with associated methods including
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import sys
import re
import numpy as np
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
# ~~> dependencies towards other pytel/modules
from utils.files import get_file_content
from utils.exceptions import TelemacException
from data_manip.extraction.parser_kenue import InS

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#
DAT_OPENH = re.compile(r'(DAMM|INSEL)')
DAT_CLOSH = re.compile(r'(DAMM|INSEL)')
DAT_FOOTER = re.compile(r'ENDE DATEI')

VAR_1INT = re.compile(r'(?P<before>[^+-]*?)(?P<number>\b(|[^a-zA-Z(,])(?:(\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')
VAR_1DBL = re.compile(r'(?P<number>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<after>.*?)\Z')
VAR_2DBL = re.compile(r'(?P<number1>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<number2>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')
VAR_3DBL = re.compile(r'(?P<number1>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<number2>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))[\s,;]*(?P<number3>[+-]?(|[^a-zA-Z(,])(?:(\d+(|\.)\d*[dDeE](\+|\-)?\d+|\d+\.\d+)(\b|[^a-zA-Z,)])))(?P<after>.*?)\Z')

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

# _____                      _______________________________________
# ____/ Toolbox for I2S/I3S /______________________________________/
#


def get_insel(file_name):
    # ~~ Get all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core = get_file_content(file_name)
    if not re.match(DAT_FOOTER, core[len(core)-1]):
        raise TelemacException(
             '\nCould not parse the following '
             'end line of the file: {}'.format(core[len(core)-1]))

    # ~~ First scan at INSEL and DAMM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is also fairly fast, so you might not need a progress bar
    core.pop(len(core)-1)
    poly = []
    vals = []
    typ = []
    npoin = 0
    iline = 0
    xyi = []
    val = []
    file_type = True
    while iline < len(core):
        proco = re.match(DAT_OPENH, core[iline])
        if proco:
            time = 0
        procc = re.match(DAT_CLOSH, core[iline])
        if procc:
            time = 1
        if proco or procc:
            iline += 1
            if xyi != []:
                poly.append(xyi)
                npoin += len(xyi)
                typ.append(time)
                vals.append(val)
                xyi = []
                val = []
        else:
            proc = re.match(VAR_3DBL, core[iline].strip())
            if proc:
                xyi.append((proc.group('number1'), proc.group('number2')))
                val.append(proc.group('number3'))
            else:
                file_type = False
                proc = re.match(VAR_2DBL, core[iline].strip())
                if proc:
                    xyi.append((proc.group('number1'), proc.group('number2')))
                    val.append('')
                else:
                    raise TelemacException(
                        '\nCould not parse the following '
                        'polyline record: {}'.format(core[iline]))
        iline += 1
    poly.append(xyi)
    npoin += len(xyi)
    typ.append(time)
    vals.append(val)

    # ~~ Second scan at INSEL and DAMM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is also fairly fast, so you might not need a progress bar
    if file_type:
        for pline in range(len(poly)):
            for iline in range(len(poly[pline])):
                x, y, val = poly[pline][iline]
                poly[pline][iline] = [float(x), float(y)]
                vals[pline][iline] = [float(val)]
            poly[pline] = np.asarray(poly[pline])
    else:
        for pline in range(len(poly)):
            for iline in range(len(poly[pline])):
                x, y = poly[pline][iline]
                poly[pline][iline] = [float(x), float(y)]
            poly[pline] = np.asarray(poly[pline])
            vals = []

    return file_type, npoin, poly, vals, typ


"""
    self.poly is a numpy object, while self.type is not.
"""


class Insel(object):

    def __init__(self, file_name):
        self.file_name = file_name
        self.file_type, self.npoin, self.poly, self.vals, self.typ = \
            get_insel(self.file_name)

    def toi2s(self, ins):
        ins.head = []
        if self.file_type:
            ins.file_type = 'i3s'
        else:
            ins.file_type = 'i2s'
        ins.npoin = self.npoin
        ins.poly = self.poly
        ins.vals = self.vals
        ins.type = self.typ
        ins.atrbut = None
        return ins

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#


__author__ = "Sebastien E. Bourban"
__date__ = "$09-Sep-2011 08:51:29$"


def main():
    """ Main function of parserJanet """

    print('\n\nInterpreting command line options\n'+'~'*72+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n
Tools for handling Janet native files in python.
    Janet and its related software (..., ) are property of Smile Consulting
        '''))
    parser.add_argument("args", nargs='*')
    options = parser.parse_args()

    if len(options.args) != 1:
        raise TelemacException(
                '\nThis program takes only one Insel type file '
                'as a time as input.\n'
                ' ... an i2s or i3s file of the same name '
                'will be created depending on the Insel content.')

    jan_file = options.args[0]
    head, _ = path.splitext(jan_file)
    insel = Insel(jan_file)

    ins = InS('')
    if insel.file_type:
        ins.file_type = 'i3s'
    else:
        ins.file_type = 'i2s'
    ken_file = head+'.'+ins.file_type
    insel.toi2s(ins)
    ins.put_content(ken_file)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
