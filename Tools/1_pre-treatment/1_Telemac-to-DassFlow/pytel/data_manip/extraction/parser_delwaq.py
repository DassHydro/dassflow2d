r"""@author TELEMAC-MASCARET Consortium

    @brief Tools for handling DELWAQ files when created by TELEMAC

"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from struct import unpack, pack
import re
import sys
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
from config import add_config_argument, update_config, CFGS
# ~~> dependencies towards other pytel/modules
from data_manip.formats.selafin import Selafin
from data_manip.formats.conlim import Conlim
from utils.progressbar import ProgressBar
from utils.files import get_file_content
from utils.exceptions import TelemacException

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#


def big2little_bot(file_name, fole_name):

    # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fle = open(file_name, 'rb')
    fole = open(fole_name, 'wb')
    print('           +> writing the surfaces-file: '+fole_name)

    # ~~ Read/Write dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ll, n_1, n_2, n_3, n_4, n_5, n_6, chk = unpack('>i6ii', fle.read(4+24+4))
    if ll != chk:
        raise TelemacException(
                '... Cannot read the first 6 INTEGER '
                'from your DELWAQ file')
    fole.write(pack('<i6ii', 4*6, n_1, n_2, n_3, n_4, n_5, n_6, 4*6))

    # ~~ Read areas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fle.seek(4, 1)
    area = np.asarray(unpack('>'+str(n_1)+'f', fle.read(4*n_1)))
    fle.seek(4, 1)
    fole.write(pack('<i', 4*n_1))
    fole.write(pack('<'+str(n_1)+'f', *(area)))
    fole.write(pack('<i', 4*n_1))

    fle.close()
    fole.close()


class DELWAQ(object):
    """
    Class to handle telemac DELWAQ output files
    """

    simplekeys = {
        "task": '',
        "geometry": '',
        "horizontal-aggregation": '',
        "minimum-vert-diffusion-used": '',
        "vertical-diffusion": '',
        "reference-time": '',
        "hydrodynamic-start-time": '',
        "hydrodynamic-stop-time": '',
        "hydrodynamic-timestep": '',
        "conversion-ref-time": '',
        "conversion-start-time": '',
        "conversion-stop-time": '',
        "conversion-timestep": '',
        "grid-cells-first-direction": '',
        "grid-cells-second-direction": '',
        "number-hydrodynamic-layers": '',
        "number-water-quality-layers": '',
        "hydrodynamic-file": '',
        "aggregation-file": '',
        "grid-indices-file": '',
        "grid-coordinates-file": '',
        "pointers-file": '',
        "lengths-file": '',
        "volumes-file": '',
        "areas-file": '',
        "flows-file": '',
        "salinity-file": '',
        "temperature-file": '',
        "vert-diffusion-file": '',
        "surfaces-file": '',
        "total-grid-file": '',
        "discharges-file": '',
        "chezy-coefficients-file": '',
        "shear-stresses-file": '',
        "walking-discharges-file": ''
                 }

    complxkeys = {
         "description": [],
         "constant-dispersion": [],
         "hydrodynamic-layers": [],
         "water-quality-layers": [],
         "discharges": []
                 }

    emptyline = re.compile(r'\s*\Z')
    comments = re.compile(r'[#]')
    var_dquot = re.compile(r'"(?P<dquot>[^"]*)"')
    var_squot = re.compile(r"'(?P<squot>[^']*)'")
    key_word = \
        re.compile(r'(?P<key>[^\s]+)\s*(?P<word>[^#]*)(?P<after>.*)\s*\Z',
                   re.I)
    key_field = re.compile(r'(?P<key>[^\s]+)\s*(?P<after>.*)\s*\Z', re.I)
    grp_word = re.compile(r'(?P<key>[^\s]*)\s*\Z', re.I)

    def __init__(self, file_name):

        # ~~> Read the steering file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not path.exists(file_name):
            raise TelemacException(
                    '... Could not file your DELWAQ file: '
                    '{}'.format(file_name))
        self.dwq_list = self.parse_dwq(get_file_content(file_name))

        # ~~> Read the geometry file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = self.dwq_list['grid-indices-file']
        if not path.exists(fle):
            raise TelemacException(
                    '...Could not find the GEO file: '
                    '{}'.format(fle))
        self.geo = Selafin(fle)
        self.npoin3 = int(self.dwq_list['grid-cells-first-direction'])
        if self.npoin3 != self.geo.npoin3:
            raise TelemacException(
                    '...In consistency in numbers with GEO file: '
                    '{} {}'.format(str(self.npoin3),
                                   str(self.geo.npoin3)))
        self.nseg3 = int(self.dwq_list['grid-cells-second-direction'])

        # ~~> Read the Conlim file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = self.dwq_list['grid-coordinates-file']
        if not path.exists(fle):
            raise TelemacException(
                    '...Could not find the Conlim file: '
                    '{}'.format(fle))
        self.conlim = Conlim(fle)

        # ~~> Time records ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.hydro0t = int(self.dwq_list['hydrodynamic-start-time'])
        self.hydroat = int(self.dwq_list['hydrodynamic-stop-time'])
        self.hydrodt = int(self.dwq_list['hydrodynamic-timestep'])
        self.hydroit = 1 + (self.hydroat-self.hydro0t) / self.hydrodt
        self.hydro00 = 0
        self.tfrom = self.hydro0t
        self.tstop = self.hydroat
        self.minvol = None

    def reset_dwq(self):
        self.hydro00 = self.hydro0t

    def minvol_dwq(self, value):
        self.minvol = float(value)

    def sample_dwq(self, tfrom, tstop):
        self.tfrom = float(tfrom)
        if self.tfrom < 0:
            self.tfrom = self.hydroat + self.tfrom + 1
        self.tstop = float(tstop)
        if self.tstop < 0:
            self.tstop = self.hydroat + self.tstop + 1
        if self.tfrom > self.tstop:
            self.tstop = self.tfrom

    def parse_dwq(self, lines):
        dwq_list = {}
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if re.match(self.comments, line):
                i += 1
                continue
            if re.match(self.emptyline, line):
                i += 1
                continue
            proc = re.match(self.grp_word, line)
            if proc:
                i += 1
                if proc.group('key').lower() in self.complxkeys:
                    end = 'end-' + proc.group('key')
                    sroc = re.match(self.key_word, lines[i].strip())
                    word = []
                    while sroc.group('key').lower() != end:
                        dval = re.match(self.var_dquot, lines[i].strip())
                        sval = re.match(self.var_squot, lines[i].strip())
                        sroc = re.match(self.key_field, lines[i].strip())
                        if dval:
                            word.append(dval.group('dquot').strip('"'))
                        elif sval:
                            word.append(sval.group('squot').strip("'"))
                        else:
                            word.append(lines[i].strip().strip("'"))
                        i += 1
                        sroc = re.match(self.key_word, lines[i].strip())
                    i += 1
                    dwq_list.update({proc.group('key').lower(): word})
                    continue
                else:
                    print('... Could not understand the following '
                          'complex key: ' + proc.group('key'))
            proc = re.match(self.key_word, line)
            if proc:
                i += 1
                if proc.group('key').lower() in self.simplekeys:
                    dval = re.match(self.var_dquot,
                                    proc.group('after').strip())
                    sval = re.match(self.var_squot,
                                    proc.group('after').strip())
                    key = proc.group('key').lower()
                    if dval:
                        dwq_list.update({key: dval.group('dquot').strip('"')})
                    elif sval:
                        dwq_list.update({key: sval.group('squot').strip("'")})
                    else:
                        dwq_list.update({key: proc.group('word').strip(" '")})
                else:
                    print('... Could not understand the following simple '
                          'key: ' + proc.group('key'))

        return dwq_list

    def big2little(self):

        fle = self.dwq_list["surfaces-file"]
        fole = path.splitext(path.basename(fle))[0] + '.qwd'
        big2little_bot(fle, fole)

        fle = self.dwq_list["lengths-file"]
        fole = path.splitext(path.basename(fle))[0] + '.qwd'
        self.big2little_nds(fle, fole)

        fle = self.dwq_list["pointers-file"]
        fole = path.splitext(path.basename(fle))[0] + '.qwd'
        self.big2little_nfx(fle, fole)

        fle = self.dwq_list["volumes-file"]
        fole = path.splitext(path.basename(fle))[0] + '.qwd'
        self.big2little_pts(fle, fole)

        fle = self.dwq_list["flows-file"]
        fole = path.splitext(path.basename(fle))[0] + '.qwd'
        self.big2little_vfx(fle, fole)

        fle = self.dwq_list["areas-file"]
        fole = path.splitext(path.basename(fle))[0] + '.qwd'
        self.big2little_are(fle, fole)

    def big2little_nds(self, file_name, fole_name):

        # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = open(file_name, 'rb')
        fole = open(fole_name, 'wb')
        print('           +> writing the lengths-file: '+fole_name)

        # ~~ Read/Write dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        nseg2 = (3*self.geo.nelem3 + self.conlim.nptfr)/2
        mbnd2 = np.count_nonzero(self.conlim.bor['lih'] != 2)
        n_3 = 2*self.geo.nplan*(nseg2 + mbnd2)

        # ~~ Read lengths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle.seek(4*4, 1)
        length = np.asarray(unpack('>'+str(n_3)+'f', fle.read(4*n_3)))
        fle.seek(4, 1)
        fole.write(pack('<iii', 4, 0, 4))
        fole.write(pack('<i', 4*n_3))
        fole.write(pack('<'+str(n_3)+'f', *(length)))
        fole.write(pack('<i', 4*n_3))

        # ~~ 3D lengths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # if n_4 != 0:
        #   fle.seek(4*n_4+8,1)

        fle.close()
        fole.close()

    def big2little_nfx(self, file_name, fole_name):

        # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = open(file_name, 'rb')
        fole = open(fole_name, 'wb')
        print('           +> writing the pointers-file: '+fole_name)

        # ~~ Read lengths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for _ in range(self.geo.nplan):
            for _ in range(self.nseg3):
                _, n_1, n_2, n_3, n_4, _ = unpack('>i4ii', fle.read(4+16+4))
                fole.write(pack('<i4ii', 4*4, n_1, n_2, n_3, n_4, 4*4))

        fle.close()
        fole.close()

    def big2little_pts(self, file_name, fole_name):

        # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = open(file_name, 'rb')
        fole = open(fole_name, 'wb')
        print('           +> writing the volumes-file: '+fole_name)

        # ~~ Read/Write dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        n_1 = self.npoin3+1
        minvol = self.minvol * np.ones(self.npoin3, dtype=np.float32)

        # ~~ Read volumes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pbar = ProgressBar(maxval=self.hydroit).start()
        for i in range(self.hydroit):
            _, itime = unpack('>ii', fle.read(4+4))
            volume = np.asarray(unpack('>'+str(self.npoin3)+'f',
                                       fle.read(4*self.npoin3)))
            volume = np.maximum(volume, minvol)
            fle.seek(4, 1)
            if itime >= self.tfrom and itime <= self.tstop:
                pbar.write('            ~> read iteration: '+str(itime), i)
                fole.write(pack('<ii', 4*n_1, itime-self.hydro00))
                fole.write(pack('<'+str(self.npoin3)+'f', *(volume)))
                fole.write(pack('<i', 4*n_1))
            else:
                pbar.write('            ~> ignore iteration: '+str(itime), i)
            pbar.update(i)
        pbar.finish()

        fle.close()
        fole.close()

    def big2little_vfx(self, file_name, fole_name):

        # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = open(file_name, 'rb')
        fole = open(fole_name, 'wb')
        print('           +> writing the flows-file: '+fole_name)

        # ~~ Read/Write dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        nseg2 = (3*self.geo.nelem3 + self.conlim.nptfr)/2
        mbnd2 = np.count_nonzero(self.conlim.bor['lih'] != 2)
        n_3 = (nseg2 + mbnd2)

        # ~~ Read volumes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pbar = ProgressBar(maxval=self.hydroit).start()
        for i in range(self.hydroit):
            _, itime = unpack('>ii', fle.read(4+4))
            vfluxes = np.asarray(unpack('>'+str(n_3)+'f', fle.read(4*n_3)))
            fle.seek(4, 1)
            if itime >= self.tfrom and itime <= self.tstop:
                pbar.write('            ~> read iteration: '+str(itime), i)
                fole.write(pack('<ii', 4*n_3, itime-self.hydro00))
                fole.write(pack('<'+str(n_3)+'f', *(vfluxes)))
                fole.write(pack('<i', 4*n_3))
            else:
                pbar.write('            ~> ignore iteration: '+str(itime), i)
            pbar.update(i)
        pbar.finish()

        fle.close()
        fole.close()

    def big2little_are(self, file_name, fole_name):

        # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fle = open(file_name, 'rb')
        fole = open(fole_name, 'wb')
        print('           +> writing the areas-file: '+fole_name)

        # ~~ Read/Write dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        nseg2 = (3*self.geo.nelem3 + self.conlim.nptfr)/2
        mbnd2 = np.count_nonzero(self.conlim.bor['lih'] != 2)
        n_3 = (nseg2 + mbnd2)

        # ~~ Read volumes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pbar = ProgressBar(maxval=self.hydroit).start()
        for i in range(self.hydroit):
            _, itime = unpack('>ii', fle.read(4+4))
            areas = np.asarray(unpack('>'+str(n_3)+'f', fle.read(4*n_3)))
            fle.seek(4, 1)
            if itime >= self.tfrom and itime <= self.tstop:
                pbar.write('            ~> read iteration: '+str(itime), i)
                fole.write(pack('<ii', 4*n_3, itime-self.hydro00))
                fole.write(pack('<'+str(n_3)+'f', *(areas)))
                fole.write(pack('<i', 4*n_3))
            else:
                pbar.write('            ~> ignore iteration: '+str(itime), i)
            pbar.update(i)
        pbar.finish()

        fle.close()
        fole.close()


# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__ = "Sebastien E. Bourban"
__date__ = "$21-Jun-2013 17:51:29$"


def main():
    """ Main function of parserDELWAQ """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+'~'*72+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n\
Tools for handling DELWAQ files when created by TELEMAC
        '''),
        usage=' (--help for help)\n---------\n       =>  '
              '%(prog)s [option] delwaq.cas \n---------')
    parser = add_config_argument(parser)
    parser.add_argument(
        "--reset", action="store_true",
        dest="areset", default=False,
        help="reset the start time to zero")
    parser.add_argument(
        "--minvol",
        dest="minvol", default='0.001',
        help="make sure there is a minimum volume")
    parser.add_argument(
        "--from",
        dest="tfrom", default="1",
        help="specify the first frame included")
    parser.add_argument(
        "--stop",
        dest="tstop", default="-1",
        help="specify the last frame included (negative from the end)")
    parser.add_argument("args", nargs="+")
    options = parser.parse_args()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for only one configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    update_config(options)
    CFGS.compute_compilation_info()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads command line arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if len(options.args) < 1:
        print('\nAt least one DELWAQ steering file name is required\n')
        parser.print_help()
        raise Exception
    file_names = options.args

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Loop over the DELWAQ files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for fle in file_names:

        # ~~> Parse DELWAQ steering file
        print('      ~> scanning your DELWAQ file: '+path.basename(fle))
        dwq = DELWAQ(fle)

        # ~~> Possible options so far
        if options.areset:
            dwq.reset_dwq()
        dwq.minvol_dwq(options.minvol)
        dwq.sample_dwq(options.tfrom, options.tstop)

        # ~~> Convert to Little Endian
        dwq.big2little()


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
