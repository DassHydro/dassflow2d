#!/usr/bin/env python
r"""@author TELEMAC-MASCARET Consortium

    @brief Handling .vtu and .pvd ASCII data files for
           particle tracking (drogues).
"""

from os import path, remove
import numpy as np
from utils.exceptions import TelemacException


def write_one_step_in_time(file_name,
                           xyz_position,
                           p_id,
                           state,
                           time_step,
                           time):
    """
    Write one time step in a .vtu file and creates its .pvd
    loader to synchronize physical time

    @param file_name (string)    Name of file (without extention)
    @param xyz_position (float)  Particule position
    @param p_id (int)            Particule id
    @param state (int)           Particule state
    @param time_step (int)       Iteration number
    @param time  (float)         Physical time

    """

    file_name_path_vtu = file_name+'_'+str(time_step)+'.vtu'
    number_of_points = (int)(len(xyz_position)/3)
    if path.exists(file_name_path_vtu):
        remove(file_name_path_vtu)
    with open(file_name_path_vtu, 'w') as fle:
        fle.write('<VTKFile type="UnstructuredGrid" version="0.1"> \n')
        fle.write('<UnstructuredGrid> \n')
        fle.write('<Piece NumberOfPoints="'
                  + str(number_of_points) +
                  '" NumberOfCells="1"> \n')
        fle.write('<PointData> \n')
        fle.write('<DataArray type="Int32" '
                  'Name="id" format="ascii"> \n')
        for i in range(number_of_points):
            fle.write(str(p_id[i])+' ')
            if (i+1) % 6 == 0 or (i+1) == number_of_points:
                fle.write('\n')
        fle.write('</DataArray> \n')
        fle.write('<DataArray type="Int32" '
                  'Name="state" format="ascii"> \n')
        for i in range(number_of_points):
            fle.write(str(state[i])+' ')
            if (i+1) % 6 == 0 or (i+1) == number_of_points:
                fle.write('\n')
        fle.write('</DataArray> \n')
        fle.write('</PointData> \n')
        fle.write('<Points> \n')
        fle.write('<DataArray type="Float32" '
                  'NumberOfComponents="3" format="ascii"> \n')
        for i in range(3*number_of_points):
            fle.write(str(xyz_position[i])+' ')
            if (i+1) % 6 == 0 or (i+1) == 3*number_of_points:
                fle.write('\n')
        fle.write('</DataArray> \n')
        fle.write('</Points> \n')
        fle.write('<Cells> \n')
        fle.write('<DataArray type="Int32" '
                  'Name="connectivity" format="ascii"> \n')
        for i in range(number_of_points):
            fle.write(str(i)+' ')
            if (i+1) % 6 == 0 or (i+1) == number_of_points:
                fle.write('\n')
        fle.write('</DataArray> \n')
        fle.write('<DataArray type="Int32" '
                  'Name="offsets" format="ascii"> \n')
        fle.write(str(number_of_points)+'\n')
        fle.write('</DataArray> \n')
        fle.write('<DataArray type="Int32" '
                  'Name="types" format="ascii"> \n')
        fle.write("1 \n")
        fle.write('</DataArray> \n')
        fle.write('</Cells> \n')
        fle.write('</Piece> \n')
        fle.write('</UnstructuredGrid> \n')
        fle.write('</VTKFile> \n')

    file_name_path_pvd = file_name+'.pvd'
    if time_step == 0:
        if path.exists(file_name_path_pvd):
            remove(file_name_path_pvd)
        with open(file_name_path_pvd, 'w') as fle:
            fle.write('<VTKFile type="Collection" version="0.1"> \n')
            fle.write('<Collection> \n')
            fle.write('<DataSet timestep="'+str(time)+'" '
                      'group="" part="0" file="'+file_name_path_vtu+'"/>\n')
            fle.write('</Collection> \n')
            fle.write('</VTKFile> \n')
    else:
        # amend with a new set
        with open(file_name_path_pvd, 'r+') as fle:
            file_text = fle.read()
            where_is = file_text.index('</Collection>')
            fle.seek(0)
            new_set = '<DataSet timestep="'+str(time)+'" '\
                      'group="" part="0" file="'+file_name_path_vtu+'"/>\n'
            fle.write(file_text[:where_is] + new_set + file_text[where_is:])


def convert_drogues_file_to_vtu(file_name_path,
                                output_file_name_path=None):
    """
    Convert drogues TELEMAC file into .vtu format (readable by ParaView)

    @param file_name_path (string) Name of the dogues file of TELEMAC (.dat)
    @param output_file_name_path (string) Path of output file
    """
    if not path.exists(file_name_path):
        raise TelemacException(
                '... could not find ASCII file for drogues: '
                '{}'.format(file_name_path))
    with open(file_name_path, 'r', encoding='utf-8') as fle:
        if output_file_name_path is None:
            output_file_name_path = file_name_path
        # Title
        line = fle.readline()
        # Labels
        line = fle.readline()
        in_3d = False
        # If Z exists, 3D format, else 2D
        if line.find("Z") != -1:
            in_3d = True
        # Meta data
        line = fle.readline().split(',')
        number_of_particules = (int)(line[2].split('=')[1])
        time = (float)(line[3].split('=')[1])
        # initialize time_step to 0, time will synchronize results
        time_step = 0
        position_xyz = np.zeros([3*number_of_particules], dtype=np.float)
        p_id = np.zeros([number_of_particules], dtype=np.int)
        state = np.zeros([number_of_particules], dtype=np.int)
        is_time_step = True
        while is_time_step:
            for i in range(number_of_particules):
                line = fle.readline().split(',')
                if in_3d:
                    position_xyz[3*i] = (float)(line[1])
                    position_xyz[3*i + 1] = (float)(line[2])
                    position_xyz[3*i + 2] = (float)(line[3])
                    p_id[i] = (int)(line[0])
                    state[i] = (int)(line[4])
                else:
                    position_xyz[3*i] = (float)(line[1])
                    position_xyz[3*i + 1] = (float)(line[2])
                    position_xyz[3*i + 2] = 0.
                    p_id[i] = (int)(line[0])
                    state[i] = (int)(line[3])
            write_one_step_in_time(output_file_name_path.split('.')[0],
                                   position_xyz,
                                   p_id,
                                   state,
                                   time_step,
                                   time)
            line = fle.readline().split(',')
            # Seeking for a new time step
            if line != '' and line[0].find("ZONE") != -1:
                time = (float)(line[3].split('=')[1])
                tmp_number_of_particules = (int)(line[2].split('=')[1])
                # not sure that resize if usefull...
                if tmp_number_of_particules != number_of_particules:
                    number_of_particules = tmp_number_of_particules
                    position_xyz = np.resize(position_xyz,
                                             [3*number_of_particules])
                    p_id = np.resize(p_id, [number_of_particules])
                    state = np.resize(state, [number_of_particules])
                time_step += 1
            else:
                is_time_step = False


def dat2vtu_parser(subparser):
    """
    Defines the arguments for argparse

    @param subparser (ArgumentParser) The parser

    @returns (ArgumentParser)  The updated parser
    """
    parser = subparser.add_parser('dat2vtu',
                                  help="Convert a drogues file\
                                        (output of telemac) "
                                  "into a VTU file (Paraview)")
    parser.add_argument(
        "drogues_file", default="",
        help="Name of the drogue file")
    # output name option
    parser.add_argument(
        dest="vtu_file", default="",
        help="Name of the output vtu file")

    return subparser
