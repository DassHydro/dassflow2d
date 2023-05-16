#!/usr/bin/env python
r"""@author TELEMAC-MASCARET Consortium

    @brief Tools for handling ASCII data files (such as CSV) in python.
"""
from __future__ import print_function
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
# ~~> dependencies towards standard python
import re
from os import path, remove
import numpy as np
from utils.exceptions import TelemacException

# _____                   __________________________________________
# ____/ Global Variables /_________________________________________/
#

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#
# same as in parser_fortran


def clean_spaces(istr):
    """
    Clean up string of spaces

    @param istr (string) The string to update

    @return (string) The updated string
    """
    istr = istr.strip()
    while '  ' in istr:
        istr = istr.replace('  ', ' ')
    while ' ,' in istr:
        istr = istr.replace(' ,', ',')
    while ', ' in istr:
        istr = istr.replace(', ', ',')
    return istr

# _____                  ___________________________________________
# ____/ Primary Classes /__________________________________________/
#


VAR_BRACKS = re.compile(r'(?P<name>.*?)(?P<unit>\([\w,*\s+-/:]*?\))')
CSV_HEADER = re.compile(r'[#]')


class CSV(object):

    def __init__(self, file_name='', delimiter=',', comments='#'):
        """
        Initialise the CSV class

        @param file_name (string) If given data is read from it
        @param delimiter (string) csv delimiter character
        @param comments (string) Comments character
        """
        self.rowheader = []
        self.rowvars = []
        self.rowunits = []
        self.colvars = ''
        self.colunits = ''
        self.colcore = None
        self.char = delimiter
        self.comments = comments
        if file_name != '':
            self.get_file_content(file_name)

    def get_columns(self, vrs):
        """
        Returns the values for columns having the name contained in vrs

        @param vrs (string) list of variables (name:type) to output ; separated

        @returns ???
        """
        # ~~> identify column subset
        subset = []
        colvar = self.colvars+' '+self.colunits
        rowvar = []
        out_var = vrs.split(';')
        for var in out_var:
            var_name = var.split(':')[0]
            for i in range(len(self.rowvars)):
                if var_name.lower() in self.rowvars[i].lower():
                    subset.append(i+1)
                    rowvar.append(self.rowvars[i]+' '+self.rowunits[i])
        return (colvar, self.colcore[0]), [('', rowvar, self.colcore[subset])]

    def add_columns(self, x, yval):
        """
        Add a new columns to the csv file

        @param x columns ???
        @param yval data ???
        """
        if self.colcore is None:
            xunit = '(-)'
            xname, x_0 = x
            proc = re.match(VAR_BRACKS, xname)
            if proc:
                xname = proc.group('name').strip()
                xunit = proc.group('unit').strip()
            self.colvars = xname
            self.colunits = xunit
            self.colcore = np.array([x_0])
        elif len(x[1]) != len(self.colcore[0]):
            raise TelemacException(
                    '... cannot aggregate columns of different supports: '
                    '{}'.format(repr(x[0])))
        u_0 = '(-)'
        ynames, y_0 = yval
        dim = len(ynames) - 1
        # ynames[0] is an extra meta data
        if dim == 1:  # /!\ This has been checked
            # ~~> This is called for cast objects, for intance
            n_0 = ynames[1]
            for i_0 in range(len(n_0)):  # each variables
                proc = re.match(VAR_BRACKS, n_0[i_0])
                if proc:
                    n_0[i_0] = proc.group('name').strip()
                    u_0 = proc.group('unit').strip()
                self.rowvars.append(n_0[i_0])
                self.rowunits.append(u_0)
                self.colcore = np.vstack((self.colcore, y_0[i_0]))
        elif dim == 2:
            # ~~> This is called for 1d:history, for instance
            n_0, n_1 = ynames[1:]
            for i_1 in range(len(n_1)):
                for i_0 in range(len(n_0)):
                    proc = re.match(VAR_BRACKS, n_0[i_0])
                    if proc:
                        n_0[i_0] = proc.group('name').strip()
                        u_0 = proc.group('unit').strip()
                    self.rowvars.append(n_0[i_0]+':'+str(n_1[i_1]))
                    self.rowunits.append(u_0)
                    self.colcore = np.vstack((self.colcore, y_0[i_0][i_1]))
        elif dim == 3:       # /!\ This has been checked
            # ~~> This is called for 1d:v-section, for instance
            n_0, n_1, n_2 = ynames[1:]
            for i_2 in range(len(n_2)):       # each plan
                for i_1 in range(len(n_1)):    # each time
                    for i_0 in range(len(n_0)):  # each variables
                        self.rowvars.append(n_0[i_0]+':'+str(n_1[i_1])+'_' +
                                            str(n_2[i_2]))
                        self.rowunits.append(u_0)
                        self.colcore = np.vstack((self.colcore,
                                                  y_0[i_0][i_1][i_2]))
        elif dim == 4:
            n_0, n_1, n_2, n_3 = ynames[1:]
            for i_3 in range(len(n_3)):
                for i_2 in range(len(n_2)):
                    for i_1 in range(len(n_1)):
                        for i_0 in range(len(n_0)):
                            self.rowvars.append(n_0[i_0]+':'+str(n_1[i_1]) +
                                                '_'+str(n_2[i_2]) +
                                                '_'+str(n_3[i_3]))
                            self.rowunits.append(u_0)
                            self.colcore = np.vstack((self.colcore,
                                                      y_0[i_0][i_1][i_2][i_3]))

    def put_file_content(self, file_name):
        """
        Write content of class into a file.

        @param file_name Name of the file that will contain the data
        """
        if path.exists(file_name):
            remove(file_name)
        with open(file_name, 'wb') as fle:
            if len(self.rowheader) > 0:
                fle.write('\n'.join(self.rowheader)+'\n')
            else:
                fle.write(b'#\n#\n')  # TODO: use header for meta data
            fle.write('{},{}\n'.format(self.colvars, ','.join(self.rowvars))
                      .encode('utf-8'))
            fle.write('{},{}\n'.format(self.colunits, ','.join(self.rowunits))
                      .encode('utf-8'))
            np.savetxt(fle, self.colcore.T, delimiter=',')

    def get_file_content(self, file_name):
        """
        Read data from an ascii file

        @param file_name (string) Name of the input file
        """
        if not path.exists(file_name):
            raise TelemacException(
                    '... could not find your ASCII file: '
                    '{}'.format(file_name))
        with open(file_name, 'r', encoding='utf-8') as fle:
            # ~~> parse header
            is_head = True
            while is_head:
                line = fle.readline()
                if line[0] == '#':
                    self.rowheader.append(line.rstrip())
                else:
                    is_head = False
            # ~~> parsing character
            if self.char not in line.rstrip():
                self.char = ' '

            # ~~> column header - variables
            vrs = clean_spaces(line).split(self.char)
            self.colvars = vrs[0]
            self.rowvars = vrs[1:]
            # ~~> column header - unit
            line = fle.readline()
            units = clean_spaces(line).split(self.char)
            self.colunits = units[0]
            self.rowunits = units[1:]
            if len(self.rowvars) != len(self.rowunits):
                raise TelemacException(
                        '... variables and units are incorrectly '
                        'numbered in {}'.format(file_name))
            # ~~> parse main values
            data = np.loadtxt(fle, comments=self.comments, delimiter=self.char)
            self.colcore = data.T
