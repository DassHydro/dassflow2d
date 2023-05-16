r"""@author TELEMAC-MASCARET Consortium

    Conversion of the i2s/i3s => shape files.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards other modules
from data_manip.extraction.parser_kenue import InS
from data_manip.formats import shapefile
# ~~> dependencies towards standard python
from os import path
# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#


class Ins2shp(InS):

    def __init__(self, file_name):
        InS.__init__(self, file_name)

    def put_content(self, file_name):
        f = shapefile.Writer(target=file_name)
        f.field('id', 'C', '40')
        i = 0                     # attributes starts at one
        for poly in self.poly:
            if self.type[i] == 1:
                f.poly([poly])
            else:
                f.line([poly])
            atrs = []
            for k in self.atrbut:
                atrs.append(self.atrbut[k][i])
            f.record(*atrs)
            i += 1
        f.close()
        return

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#


__author__ = "Sebastien E. Bourban"
__date__ = "$13-Jan-2012 08:51:29$"


def kenue2shp(ins_files):
    """
    """
    for ins_file in ins_files:

        ins_file = path.realpath(ins_file)
        print('\n\nSmoothing ' + path.basename(ins_file) + ' within ' +
              path.dirname(ins_file) + '\n'+'~'*72+'\n')
        ins = Ins2shp(ins_file)

        ins.put_content(ins_file)


def kenue2shp_parser(subparser):
    """
    Adding options for kenue2srf

    @param subparser (ArgumentParser) The argument parser to update

    @return (ArgumentParser) the updated argument parser
    """

    parser = subparser.add_parser('kenue2shp',
                                  help='Convert a Kenue file into a shapefile')
    parser.add_argument("args", nargs="+")

    return subparser
