r"""@author Sebastien E. Bourban

    @brief
        Tools for handling MSH files when created by the mesh generator GMSH

    @details
        Contains read/write functions for binary and asci MSH files
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
from struct import unpack
import re
import sys
from os import path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
from matplotlib.tri import Triangulation
# ~~> dependencies towards other pytel/modules
from data_manip.formats.selafin import Selafin
from data_manip.extraction.parser_kenue import InS
from utils.files import put_file_content
from utils.exceptions import TelemacException

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#


class GEO(InS):

    def __init__(self, file_name):

        # ~~> possibly empty i2s
        InS.__init__(self, file_name)

        # TODO: You may need to account for the fact that more than one polygon
        #       is an outside domain boundary
        # ~~> sort out by the area to build the order map
        self.sort_by_areas()
        # ~~> make the first / larger polygon anti-clockwise
        self.make_anti_clockwise(select=[0])
        # ~~> make the other polygons clockwise
        self.make_clockwise(select=range(len(self.poly))[1:])

        # ~~> initialisasing counters
        self.ipoin = [0]
        self.iline = [0]
        self.iloop = [0]
        self.isurf = [0]

    def write_polygon(self, poly):
        # TODO: values / resolution
        geo = []
        ipoin = self.ipoin[-1]
        for i, x_y in zip(range(len(poly)), poly):
            geo.append('Point('+str(i+ipoin+1)+') = { ' +
                       str(x_y[0])+','+str(x_y[1])+',0,2000 };')
        # Lines
        iline = self.iline[-1]
        for i in range(len(poly))[:-1]:
            geo.append('Line('+str(i+iline+1)+') = { ' +
                       str(i+iline+1)+','+str(i+iline+2)+' };')
        geo.append('Line('+str(len(poly)+iline)+') = { ' +
                   str(len(poly)+iline)+','+str(iline+1)+' };')
        # Line Loop
        iloop = self.iloop[-1]
        geo.append('Line Loop('+str(len(poly)+iloop+1)+') = {' +
                   ','.join([str(i+iloop+1) for i in range(len(poly))])+' };')
        # next set of entities
        self.ipoin.append(ipoin+len(poly))
        self.iline.append(iline+len(poly)+1)
        self.iloop.append(iloop+len(poly)+1)
        return geo

    def put_content(self, file_name, head=None):

        geo = []
        # assuming the first polygon is the outside domain
        i_p = 0
        geo.extend(self.write_polygon(self.poly[i_p]))
        # add the other polygons
        for i_p in range(1, self.npoly):
            geo.extend(self.write_polygon(self.poly[i_p]))

        # make up surface with wholes
        psurf = self.isurf[-1]
        geo.append('Plane Surface('+str(psurf+1)+') = {' +
                   ','.join([str(i) for i in self.iloop[1:]])+'};')
        self.isurf.append(psurf+1)

        # write up
        put_file_content(file_name, geo)


class MSH(Selafin):

    mshkeys = {"MeshFormat": '',
               "Nodes": '',
               "Elements": [],
               "PhysicalName": '',
               "Periodic": '',
               "NodeData": '',
               "ElementData": '',
               "ElementNodeData": '',
               "InterpolationScheme": ''}

    frst_keys = re.compile(r'[$](?P<key>[^\s]+)\s*\Z', re.I)
    last_keys = re.compile(r'[$]End(?P<key>[^\s]+)\s*\Z', re.I)

    def __init__(self, file_name):

        # ~~> empty Selafin
        Selafin.__init__(self, '')
        self.datetime = []

        # ~~> variables
        self.title = ''
        self.nbv1 = 1
        self.nvar = self.nbv1
        self.varindex = range(self.nvar)
        self.varnames = ['BOTTOM          ']
        self.varunits = ['M               ']

        # ~~ Openning files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.fle = {}
        self.fle.update({'name': file_name})
        # "<" means little-endian, ">" means big-endian
        self.fle.update({'endian': ">"})
        self.fle.update({'integer': ('i', 4)})  # 'i' size 4
        self.fle.update({'float': ('f', 4)})  # 'f' size 4, 'd' = size 8
        self.fle.update({'hook': open(file_name, 'r')})
        fle = iter(self.fle['hook'])

        # ~~ Read/Write dimensions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Note:
        #    The section MeshFormat is mandatory
        line = fle.next()
        proc = re.match(self.frst_keys, line)
        if proc:
            if proc.group('key') != "MeshFormat":
                raise TelemacException(
                    '... Could not recognise your MSH file format. '
                    'Missing MeshFormat key.')
            line = fle.next().split()
            if line[0] != "2.2":
                raise TelemacException(
                    '... Could not read your MSH file format. '
                    'Only the version 2.2 is allowed.')
            file_type = int(line[1])
            if file_type == 1:
                print('... I have never done this before. Do check it works')
                line = fle.next()
                _, _, _ = unpack('>i', line.read(4+4+4))
            float_size = int(line[2])
            if float_size == 8:
                self.fle['float'] = ('d', 8)
        line = fle.next()
        proc = re.match(self.last_keys, line)
        if proc:
            if proc.group('key') != "MeshFormat":
                raise TelemacException(
                        '... Could not complete reading the header of you MSH '
                        'file format. Missing EndMeshFormat key.')

        # ~~ Loop on sections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        while True:
            try:
                line = fle.next()
            except StopIteration:
                break
            proc = re.match(self.frst_keys, line)
            if not proc:
                raise TelemacException(
                    '... Was expecting a new Section starter. '
                    'Found this instead: {}'.format(line))
            key = proc.group('key')

        # ~~ Section Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if key == "Nodes":
                print('     +> mesh x,y,z')
                npoin = int(fle.next())
                if self.fle['float'][0] == 'd':
                    meshx = np.zeros(npoin, dtype=np.float64)
                    meshy = np.zeros(npoin, dtype=np.float64)
                    meshz = np.zeros(npoin, dtype=np.float64)
                else:
                    meshx = np.zeros(npoin, dtype=np.float)
                    meshy = np.zeros(npoin, dtype=np.float)
                    meshz = np.zeros(npoin, dtype=np.float)
                # map_nodes = []
                for i in range(npoin):
                    line = fle.next().split()
                    # map_nodes.append(int(line[0]))
                    meshx[i] = np.float(line[1])
                    meshy[i] = np.float(line[2])
                    meshz[i] = np.float(line[3])
                # TODO: renumbering nodes according to map_nodes ?
                # map_nodes = np.asarray(map_nodes)
                self.npoin2 = npoin
                self.meshx = meshx
                self.meshy = meshy
                self.meshz = meshz

                line = fle.next()

        # ~~ Section Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            elif proc.group('key') == "Elements":
                print('     +> renumbered connectivity')
                nelem = int(fle.next())
                ikle2 = - np.ones((nelem, 3), dtype=np.int)
                for i in range(nelem):
                    line = fle.next().split()
                    if int(line[1]) != 2:
                        continue
                    expr = line[int(line[2])+3:]
                    ikle2[i] = [np.int(expr[0]), np.int(expr[1]),
                                np.int(expr[2])]

                self.ikle2 = ikle2[np.not_equal(*(np.sort(ikle2).T[0::2]))] - 1
                self.nelem2 = len(self.ikle2)

                line = fle.next()
                # TODO: fitting the unique node numbers with map_nodes ?

        # ~~ Unnecessary section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            else:
                while True:
                    line = fle.next()
                    if re.match(self.last_keys, line):
                        break

        proc = re.match(self.last_keys, line)
        if proc:
            if proc.group('key') != key:
                raise TelemacException(
                        '... Could not complete reading the header of your '
                        'MSH file format. Missing {} end key.'.format(key))

        # ~~> sizes
        print('     +> sizes')
        self.ndp3 = 3
        self.ndp2 = 3
        self.nplan = 1
        self.nelem3 = self.nelem2
        self.npoin3 = self.npoin2
        self.ikle3 = self.ikle2
        self.iparam = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

        print('     +> boundaries')
        # ~~> establish neighborhood
        _ = Triangulation(self.meshx, self.meshy, self.ikle3)\
            .get_cpp_triangulation().get_neighbors()
        # ~~> build the enssemble of boundary segments
        # ~~> define ipobO from an arbitrary start point
        self.ipob3 = np.ones(self.npoin3, dtype=np.int)
        self.ipob2 = self.ipob3

    def put_content(self, file_name, showbar=True):

        # ~~> new SELAFIN writer
        self.fole = {}
        self.fole.update({'hook': open(file_name, 'wb')})
        self.fole.update({'name': file_name})
        self.fole.update({'endian': ">"})     # big endian
        self.fole.update({'float': ('f', 4)})  # single precision

        print('     +> Write SELAFIN header')
        self.append_header_slf()

        print('     +> Write SELAFIN core')
        self.append_core_time_slf(0.0)
        self.append_core_vars_slf([self.meshz])
        self.fole['hook'].close()

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#


__author__ = "Sebastien E. Bourban"
__date__ = "$11-Nov-2015 17:51:29$"


def main():

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nLoading Options and Configurations\n'+'~'*72+'\n')
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=('''\n\
Tools for handling MSH files when created by the mesh generator GMSH
        '''))
    parser.add_argument("args", nargs='*')
    options = parser.parse_args()
    if len(options.args) < 1:
        raise TelemacException(
             '\nThe name of and action is required, together with '
             'associated arguments\n')

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Reads code name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    code_name = options.args[0]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # ~~~~ Case of MSH to SELAFIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if code_name == 'msh2slf':

        # ~~ Reads command line arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if len(options.args) != 2:
            print('.. One MSH file name is required\n')
            parser.print_help()
            raise TelemacException
        file_name = options.args[1]
        if not path.exists(file_name):
            raise TelemacException(
                    '... Could not file your MSH file: '
                    '{}'.format(file_name))

        # ~~ Parse the MSH file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('      ~> scanning your MSH file: '+path.basename(file_name))
        msh = MSH(file_name)
        head, _ = path.splitext(file_name)
        # ~~ Convert to SELAFIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('      ~> converting it to a SELAFIN: '+path.basename(file_name))
        msh.put_content(head+'.slf')

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Case of i2s/i3s to GEO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elif code_name == 'ins2geo':

        # ~~ Reads command line arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if len(options.args) != 2:
            print('... One i2s/i3s file name is required\n')
            parser.print_help()
            raise TelemacException
        file_name = options.args[1]
        if not path.exists(file_name):
            raise TelemacException(
                    '... Could not file your i2s/i3s file: '
                    '{}'.format(file_name))

        # ~~ Parse the i2s/i3s file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('      ~> scanning your MSH file: {}'
              .format(path.basename(file_name)))
        geo = GEO(file_name)
        head, _ = path.splitext(file_name)
        # ~~ Convert to GEO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print('      ~> converting it to a SELAFIN: {}'
              .format(path.basename(file_name)))
        geo.put_content(head+'.geo')

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Case of UNKNOWN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else:
        raise TelemacException(
                '\nDo not know what to do with this code name: '
                '{}'.format(code_name))


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Jenkins' success message ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('\n\nMy work is done\n\n')

    sys.exit(0)


if __name__ == "__main__":
    main()
